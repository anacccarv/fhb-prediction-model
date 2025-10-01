
p_hat <-  predict(meta_model, type = "response")
y  <- df_predictors$epidemic
sev_epid <- data %>% 
  group_by(study) %>% 
  summarise(index = mean(index)) %>% 
  mutate(index = index/100) %>% 
  filter(index > 0.1) %>% 
  pull(index)

# ============================================================
# NMB vs threshold (p_t) — Global model (Duffeck et al. 2020)
# With risk-conditional severity, NMB per treated unit, and PPV vs C/B
# All monetary values in USD/ha (1 BRL = 0.1878 USD)
# ============================================================

# Conversion rate
brl_to_usd <- 0.1878

# Environment prerequisites:
# y (0/1 epidemic), p_hat (predicted probability), length(y)==length(p_hat)
# One of the severity options:
#   A) sev ........ vector aligned with unit (fraction 0–1; NA when no epidemic)
#   B) sev_epid ... vector of severities (fraction 0–1) only from years/locations WITH epidemic
stopifnot(exists("y"), exists("p_hat"), length(y) == length(p_hat))
have_aligned_sev <- exists("sev") && length(sev) == length(y)

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(purrr); library(tidyr); library(patchwork)
  if (have_aligned_sev) library(mgcv)
})

# ----------------------------
# 1) Problem parameters
# ----------------------------
# Duffeck et al. 2020 (global FHB model; spring wheat/Brazil)
slope_kg <- 49.1      # kg/ha lost per 1 point of FHB index

# Efficacy (more conservative; adjust if you have a meta-analysis)
rE <- function(S) rbeta(S, 40, 60)   # mean ≈ 0.40; 95% CI ~ 0.25–0.56

# Total cost (product + spraying + operation), converted to USD
rC <- function(S) pmax(rnorm(S, mean = 150, sd = 20), 0) * brl_to_usd

# Wheat price (USD/kg) — converted
P_usd_per_kg <- 1.13 * brl_to_usd
preco_saca_usd <- 60 * P_usd_per_kg   # 60 kg sack

# Threshold grid and PSA size
ths <- seq(0.05, 0.60, by = 0.01)
S   <- 10000
set.seed(123)

# ----------------------------
# 2) Risk-conditional severity  E[s_pts | p]
# ----------------------------
h_bw <- 0.06

if (have_aligned_sev) {
  base_idx <- which(y == 1 & is.finite(p_hat) & is.finite(sev))
  p_ep  <- p_hat[base_idx]
  s_ep  <- sev[base_idx]              # fraction 0–1
  draw_s_pts <- function(pt, S){
    w <- exp(- (p_ep - pt)^2 / (2*h_bw^2)) 
    if (sum(w) == 0 || all(!is.finite(w))) w <- rep(1, length(p_ep))
    w <- pmax(w, 1e-9)
    s_frac <- sample(s_ep, size = S, replace = TRUE, prob = w)
    100 * s_frac  # return in POINTS
  }
} else {
  stopifnot(exists("sev_epid"))
  s_base <- sev_epid[is.finite(sev_epid)]
  s_base <- s_base[s_base >= 0 & s_base <= 1]
  k_prec <- 25
  q <- seq(0, 1, by = 0.1)
  mu_q <- stats::quantile(s_base, probs = q, na.rm = TRUE)
  fit_mu <- stats::splinefun(x = q, y = mu_q, method = "monoH.FC")
  draw_s_pts <- function(pt, S){
    r <- mean(p_hat[y==1] <= pt, na.rm = TRUE)  
    m  <- min(max(fit_mu(r), 0), 1)            
    alpha <- max(m * k_prec, 1e-3); beta <- max((1 - m) * k_prec, 1e-3)
    100 * rbeta(S, alpha, beta)                 
  }
}
k_prec <- 20

# ----------------------------
# 3) PSA per threshold (population NMB and per treated unit)
# ----------------------------
compute_psa <- function(y, p_hat, ths, slope_kg, P, rE, rC, S = 10000){
  N <- length(y)
  map_dfr(ths, function(pt){
    pred <- p_hat >= pt
    TP <- sum(pred & y == 1)
    Tt <- sum(pred)
    PPV <- if (Tt > 0) TP / Tt else 0
    
    # Sampling
    E  <- rE(S)
    C  <- rC(S)
    s_pts <- draw_s_pts(pt, S)
    
    # Benefit per epidemic (USD/ha): slope_kg * E * s_pts * P
    B  <- slope_kg * E * s_pts * P
    
    # NMB per decision unit (population) and per treated unit
    NMB_pop  <- (TP / N) * B - (Tt / N) * C
    NMB_trat <- PPV * B - C
    
    req_PPV  <- pmin(C / pmax(B, 1e-9), 1)
    
    tibble(
      pt,
      NMB_mean = mean(NMB_pop),
      NMB_lo   = unname(quantile(NMB_pop, .025)),
      NMB_hi   = unname(quantile(NMB_pop, .975)),
      Pr_pos   = mean(NMB_pop > 0),
      Pr_ge_1saca = mean(NMB_pop >= preco_saca_usd),
      NMB_trat_mean = mean(NMB_trat),
      Pr_trat_pos   = mean(NMB_trat > 0),
      PPV = PPV,
      PPV_req_med = median(req_PPV),
      PPV_req_lo  = quantile(req_PPV, .10),
      PPV_req_hi  = quantile(req_PPV, .90)
    )
  })
}

res <- compute_psa(y, p_hat, ths, slope_kg, P_usd_per_kg, rE, rC, S)

# ----------------------------
# 4) Plots
# ----------------------------
pA <- ggplot(res, aes(pt, NMB_mean)) +
  geom_ribbon(aes(ymin = NMB_lo, ymax = NMB_hi), alpha = 0.15) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Risk threshold (p_t)", y = "Mean NMB (USD/ha per decision unit)",
       title = "(A) Mean NMB and 95% CI") +
  theme_minimal()

pB <- ggplot(res, aes(pt, Pr_pos)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_hline(yintercept = 0.8, linetype = 3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
  labs(x = "Risk threshold (p_t)", y = "Probability of NMB > 0",
       title = "(B) Cost-effectiveness acceptability curve (Pr[NMB>0])") +
  theme_minimal()

pC <- ggplot(res, aes(pt, Pr_ge_1saca)) +
  geom_line(size = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
  labs(x = "Risk threshold (p_t)", y = "Probability of NMB ≥ 1 sack (60 kg)",
       title = "(C) Minimum gain target (≥ 1 60kg bag)") +
  theme_minimal()

pD1 <- ggplot(res, aes(pt, NMB_trat_mean)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Risk threshold (p_t)", y = "NMB per treated unit (USD/ha)",
       title = "(D1) NMB per treated unit = PPV·B − C") +
  theme_minimal()

pD2 <- ggplot(res, aes(pt, PPV)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = PPV_req_lo, ymax = PPV_req_hi), alpha = 0.15, fill = "grey70") +
  geom_line(aes(y = PPV_req_med), linetype = 2) +
  scale_y_continuous(limits = c(0,1), labels = scales::percent) +
  labs(x = "Risk threshold (p_t)", y = "PPV (observed) vs PPV* (C/B)",
       title = "(D2) Required PPV (10–90% range) and observed") +
  theme_minimal()

(pA | pB) /  (pD1 | pD2)
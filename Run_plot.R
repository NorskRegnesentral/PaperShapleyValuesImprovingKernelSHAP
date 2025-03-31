# This file assumes that Run_sampling.R and Run_SV.R has been run for the given number of rho values and repetitions
library(data.table)
library(ggplot2)

# The parameters of the experiment
m = 10
n_train = 1000
n_test = 1000
rhos = c(0, 0.2, 0.5, 0.9)
rho_equi = TRUE
betas = c(2, 10, 0.25, -3, -1, 1.5, -0.5, 10, 1.25, 1.5, -2)

# Where rds objects are saved
folder_save = "./rds"

# The number of repetitions. Change to seq(500) for 500 repetitions
repetitions = 1

# Load the MAE data
res = data.table::rbindlist(
  lapply(rhos, function(rho) {
    dt_now_tmp = data.table::rbindlist(
      lapply(repetitions, function(repetition) {
        # Get the true value
        file_name = file.path(folder_save,
                              paste0(paste("M", m, "n_train", n_train, "n_test", n_test,  "rho", rho, "equi", rho_equi,
                                           "betas", paste(as.character(betas), collapse = "_"),
                                           "MAE_repetition", repetition, sep = "_"),
                                     ".rds"))
        if (file.exists(file_name)) {

          res_tmp = readRDS(file_name)
        } else {
          res_tmp = NULL
        }
        return(res_tmp)
      }
      )
    )
  }
  )
)

# Compute the mean, lower, and upper MAE
res_MAE = res[, .(MAE_mean = mean(MAE),
                  MAE_lower = quantile(MAE, 0.025),
                  MAE_upper = quantile(MAE, 0.975)),
              by = c("Rho", "N_S", "Strategy")]

# Make a plot of the MAE with the emprical confidence bands
fig_MAE =
  ggplot(res_MAE, aes(x = N_S, y = MAE_mean, col = Strategy, fill = Strategy)) +
  facet_wrap( . ~ Rho, labeller = label_bquote(cols = rho ==.(Rho)), scales = "free_y") +
  geom_ribbon(aes(ymin = MAE_lower, ymax = MAE_upper), alpha = 0.4, linewidth = 0.0) +
  geom_line(linewidth = 0.65) +
  scale_x_continuous(labels = scales::label_number()) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  guides(col = guide_legend(nrow = 2, byrow = TRUE),
         fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(color = "Strategy:", fill = "Strategy:", linetype = "Strategy:",
       x = expression(N[S]),
       y = bquote(bar(MAE)[500]*"("*bold(phi)*", "*bold(phi)[italic(D)]*")")) +
  theme_minimal() +
  theme(legend.position = 'bottom',
        strip.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.6)),
        axis.text = element_text(size = rel(1.5)))

fig_MAE

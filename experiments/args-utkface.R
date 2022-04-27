# Arguments used to run utkface experiments
# Andrea Goetschi
# April 2022

# Arguments were adjusted for each experiment

fname <- "utkface_ci_lossrps_wsyes_augno"
mod <- "ci"
# fml <- age_group ~ gender + x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10
fml <- age_group ~ 1
# loss <- "nll"
loss <- "rps"
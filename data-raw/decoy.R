## code to prepare `decoy` dataset goes here

decoy <- read.csv("http://datacolada.org/appendix/74/Study 3 - Decoy Effect.csv")

decoy <- decoy %>%
  tidyr::pivot_longer(cols = Day1:Day40,
                      names_to = "day",
                      values_to = "weight") %>%
  dplyr::select(-Group..1.experimental.condition..2..control.condition.,
                -day) %>%
  dplyr::rename(subject = Subject,
                workroom = WorkroomNo.)


usethis::use_data(decoy, overwrite = TRUE)

library(dplyr)
library(ggplot2)

setwd("/home/ecs-assist-user/loadtest/app_v0.2/all.v2")

df <- shinyloadtest::load_runs(
  `1U` = "/home/ecs-assist-user/loadtest/app_v0.2/all/run1",
  `10U` = "/home/ecs-assist-user/loadtest/app_v0.2/all/run10",
  `20U` = "/home/ecs-assist-user/loadtest/app_v0.2/all/run20",
  `30U` = "/home/ecs-assist-user/loadtest/app_v0.2/all/run30",
  `40U` = "/home/ecs-assist-user/loadtest/app_v0.2/all/run40",
  `50U` = "/home/ecs-assist-user/loadtest/app_v0.2/all/run50",
  `60U` = "/home/ecs-assist-user/loadtest/app_v0.2/all/run60"
)
shinyloadtest::shinyloadtest_report(df, "shinyloadtest_report.html")
colnames(df)

df$run %>% table

df$user_id %>% table

df$session_id %>% table

df$iteration %>% table

df$input_line_number %>% table

df$event %>% table

df$start
df$end
df$time
df$concurrency

df$maintenance %>% table
df$label %>% table
df$label %>% class()
df$json

events = df %>% filter(maintenance == T) %>% group_by(label) %>% summarise(avg=mean(time)) %>% arrange(desc(avg)) %>% head(8) %>% select(label) %>% dplyr::pull()
events

result = df %>% filter(df$label %in% events & maintenance == T) %>% filter(run %in% c("1U","10U","20U","30U","40U","50U"))
ggplot(result, aes(x=run, y = time, fill = run)) + 
  facet_wrap(~ label, ncol = 4) + 
  geom_boxplot() + 
  labs(x = "", y = "Time (sec)") + 
  theme(
    legend.position = "none"
  )
ggsave('event_time.pdf', width = 9, height = 4.5)


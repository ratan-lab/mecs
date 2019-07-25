#!/usr/bin/env Rscript 

library(tidyverse)

muts <- read_tsv("mutations.tsv", col_names=F)
found <- read_tsv("variants.ann.txt", col_names=F) %>% 
         select(X2,X3) %>%
         mutate(X3 = X3-1) %>% 
         distinct() %>% mutate(found = 1)

pdf("assessment.pdf")
muts %>% left_join(found, by=c("X1"="X2", "X2"="X3")) %>%
    mutate_all(funs(replace_na(., 0))) %>%
    mutate(found = factor(found)) %>% 
    ggplot(aes(X6,X3,color=found)) + geom_point()
dev.off()

---
layout: post
title: Creating a duty schedule in R
date: 2019-08-06
categories: Programming
tags: R
image: /assets/posts/2019-08-06-duty-schedule-in-r/duty-schedule-in-r_thumbnail.jpeg
_focus_key_word: duty schedule in R
excerpt: "By creating a duty schedule in R you avoid bias and save time. You just provide a list of people and get a PDF of the duty schedule."
---

As a person who possesses some programming skills, I try to automate everything I can. Recently, I became responsible for creating a kitchen duty schedule at work. So, I wrote an **R script** that takes a list of people as input and outputs a PDF with the schedule and I would like to share it with you.

## Schedule requirements

The duty assumes one person cleans the **kitchen** for a week and another person makes **[fika](https://en.wikipedia.org/wiki/Coffee_culture#Sweden){:target="_blank"}** on that week. It is also essential to take into account that **the same person should not be responsible for both** kitchen and fika during the same year. The frequency of being in the schedule list should also be **fairly distributed among people**.

## Generate a schedule table

First, you need to **load the list** of people, and extract the names:

```r
d <- read.table("people-list.csv", sep = "\t", header = T, stringsAsFactors = F)
names <- d$Name
```

Then, **randomly pick** a few people (in my case it was 9) who will be assigned to the kitchen duty:

```r
Kitchen <- sample(names, 9)
```

Do the same to assign the fika duty but make sure people from the kitchen duty list are **excluded**:
```r
Fika <- sample(names[!(names %in% Kitchen)], 9)
```

After the people lists are created, generate the **start and end dates** as well as **week numbers** for these lists:

```r
start <- seq(as.Date("19/08/19", format = "%d/%m/%y"), by = "week", length.out = 9)
end <- seq(as.Date("23/08/19", format = "%d/%m/%y"), by = "week", length.out = 9)
week <- strftime(start, format = "%V")
```

In the end, merge these list into a **table**:
```r
dd <- data.frame(Kitchen, Fika, start, end, week)
write.table(dd, "kitchen-schedule_week34-42.csv", sep = "\t", row.names = F)
```
Everything seems to be done. One could just load this table into a spreadsheet editor, format it to a nice look and print. But why waste time on this manual work if you can automate this step too.

## Plot a table in R

Instead of manually formatting the obtained table in a spreadsheet editor, you can add a few more lines R code and get a **print-ready table**:

```r
library(gridExtra)
library(grid)

pdf("kitchen-schedule_week34-42.pdf", width=11.69, height=8.27)
g <- tableGrob(dd, rows = NULL, theme = ttheme_default(base_size = 16,
               padding = unit(c(20, 12), "mm")))
grid.newpage()
grid.draw(g)
dev.off()
```

In the end, you will obtain a PDF page of **A4 size** with this kind of table:

![A schedule table generated in R](/assets/posts/2019-08-06-duty-schedule-in-r/schedule-table-in-r.jpeg)

## Generating new schedule tables

Next time you generate a schedule table, you just need to **exclude** the people who were assigned some **duties before**:

```r
d <- read.table("people-list.csv", sep = "\t", header = T, stringsAsFactors = F)
toexlcude <- read.table('kitchen-schedule_week34-42.csv',
                        header = T, sep = "\t", stringsAsFactors = F)

names <- d$Name[!(d$Name %in% c(toexlcude$Kitchen, toexlcude$Fika))]
```

The rest of the code is the same as above. If you have several duty lists with the names you need to exclude, just merge them before applying the exclusion.


## Full Code

All the code put together:

```r
library(gridExtra)
library(grid)

d <- read.table("people-list.csv", sep = "\t", header = T, stringsAsFactors = F)

if(file.exists('previous_kitchen-schedule.csv')){
  toexlcude <- read.table('previous_kitchen-schedule.csv',
                          header = T, sep = "\t", stringsAsFactors = F)
  names <- d$Name[!(d$Name %in% c(toexlcude$Kitchen, toexlcude$Fika))]
}else{
  names <- d$Name
}

Kitchen <- sample(names, 9)
Fika <- sample(names[!(names %in% Kitchen)], 9)

start <- seq(as.Date("21/10/19", format = "%d/%m/%y"), by = "week", length.out = 9)
end <- seq(as.Date("25/10/19", format = "%d/%m/%y"), by = "week", length.out = 9)
week <- strftime(start, format = "%V")

dd <- data.frame(Kitchen, Fika, start, end, week)
write.table(dd, "kitchen-schedule_week43-51.csv", sep = "\t", row.names = F)
sample(d$Name, 1)

pdf("kitchen-schedule_week43-51.pdf", width=11.69, height=8.27)
g <- tableGrob(dd, rows = NULL, theme = ttheme_default(base_size = 16,
               padding = unit(c(20, 12), "mm")))
grid.newpage()
grid.draw(g)
dev.off()
```

## Final thoughts

If you will ever be asked to **volunteer for creating duty schedules**, do not hesitate to agree. It will cost you so little. Just modify this script for your needs and generate a duty schedule in R with one click.

*If you have any questions or suggestions, feel free to [email me](mailto:dmytro.kryvokhyzha@evobio.eu)*.

df<-read.csv("2_seminaro_užduotys.csv", header=TRUE, stringsAsFactors = FALSE)

#1
colnames(df)
#2
head(df, 4)
#3
tail(df, 4)
#4
nrow(df)
ncol(df
     )
summary(df)
#5
df$Ozone[47]
#6
df$Wind[50:59]
#7
sum(is.na(df$Ozone))
#8
mean(df$Ozone, na.rm=TRUE)
#9
mean(df$Solar.R[df$Ozone>31&df$Temp>90], na.rm = T)
#10
mean(df$Temp[df$Month==6], na.rm = T)
max(df$Temp[df$Month==5], na.rm = T)

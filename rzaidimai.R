#R turi 5 bazinius objektu tipus
#1. numeric 1,4 -5,3
#2. integr 1L 2L (sveikas skaicius)
#3. complex 1i 2i
#4. character 'vilnius' 'amzius'
#5. logical TRUE/FALSE
#su komanda typeof() galime pasitikrinti klase

#c()funckija leidzia sukurrti objektu vektorius
# x <- c(0.2 , 0.6 ) # numeric class
# x <- c(T, F) #logical class
# x <- c("a ", "b" ,"c") #character class
# x <- 1:5 #integer class
# x <- c(1+0i , 2+4i) #complex class

#norint sukurti tuscia vektoriu
# x <- vector(mode = "numeric", length = 5)
# bus sukurtas vektorius 0 0 0 0 0

#jei jungiami skirtingu klasiu objektai
#R priskiria bendriausia klase visiems vektoriuje objektams
#TRUE =1, FALSE=0
#procesas vadinamas coersion(vektoriu sujungimas)
#pvz x <- c("a", T, F) #character class

#norint priskirti tam tikrai klasei
#rasome pvz as.numeric(x), as.factor(x)
#jei bus nelogiski priskyrimai ismes NA

# list() talpina ivairiu klasiu objektus
# x<- list(1.2,TRUE,2L)

#norint sukurti matrica rasom x<- matrix(1:9,nrow=3, ncol=3)
# x<- matrix(1:9,nrow=3, ncol=3, byrow=TRUE) byrow- skaiciai eina eilutem
# dim(x) <- c(4,3) - sukuriam vektoriui dimensijas
# cbind - is skirtingu vektoriu sujungia matricas stulepliais
# rbind - vektorius raso eilute

# x<- factor(c("taip, "ne" , "taip"), levels=c("taip","ne"))
#table(x) - suskaiciuoja kiek yra vienu ir kitu pvz( girtas 3 , blaivas 2)

# is.na() testuoja ar egz NA, is.nan - egz ar yra NaN

# x<- data.frame(Fname=C("ana","maria","John), Grades=c(9,10,7))
# nrow(x) suskaiciuoja eiles, ncol(x) suskaiciuoja column, rownames ir colnames paraso pavadinimus
# priskiriant rasom colnames(x) <- "kazkas" rownames(x)<- c("alpha", "beta", "gama")
# keiciant pavadinimus matricom rasom dimnames(x)<- list(c("a", "b"), c("c","d"))

# x<- array(c(1:6),dim = c(3,3,2)) 2 gale parodo kad bus 2 masyvai

# read.table() - nuskaito faila header=TURE - kad yra stulp pavadinimai

#norint paziureti kurioj direktorijoj esi getwd() nustatyti - setwd()


#pagrindiniai operatoriai leidziantys pasirinkti dali objektu:
# [x] - visad duoda objekta tos pacios klases, galima pasirinkti daugiau nei viena elementa
# [[x]] - vieno elemento is list arba data frame pasirinkimui
# $ ledizia pasirinkti pagal pavadinimus (pagal col.names)

#norint isvalyti NA galime rasyti:
# x<- c(1,2,3,4,NA)
#isimtina <- is.na(x)
# x[!isimtina]

# x[complete.cases(x)] - parodo vektoriu be NA
# na.omit(airquality)[1:6,] - pasalina na is datasetu

#valdymo sturkturos leidzia valdyti programu veikima, nuo tam tikru aplinkybiu:
#if,else: testuoja tam tikra aplinkybe
#for: vykdo programa tam tikra iteraciju skaiciu
#while: vykdo programa kol egzistuoja tam tikros aplinkybes
#repeat: vykdo nesibaigiancia interacija
#break: nutraukia interacijos procesa
#next: persoka 1 interacija
#return: nutraukia funkcija

#if 
x<-5
if(x>3){
  y<-5
} else{
  y<-0
}

# for loop dazniausiai interuoti tam tikriems vektoriams zinant kiek kartu interacija turi trukti
for(i in 1:6){
  print(i)
}

for (values in 1:5){
  print(values)
}

for(i in 1:6) print(i)

x<- c("a", "b","c","d","e","e","f")
for (i in 1:6) {
  print(x[i])
}

#while testuoja aplinkybes. jeigu ok atlieka veiksma pavaigus vel testuoja ir t.t
# while loopai gali testis neriboda skaiciu iteraciju tad atsarigia

while(count<5000){
  print(count)
  count <- count+1
}
#cia neaisku kada baigsis interacija
z<- 5
while(z>=3 && z<=10){
  print(z)
  #rbinom(n,size,prob)
  coin <- rbinom(1,1,0.5)
  if(coin==1){
    z<-z+1
  }else{
      z<- z-1
    }
}

#repeat incializuoja begalines trukmes loop
#vienintelis budas tai sustabdyti yra su break
#pavojinga funkcija
x0<- 1
tol<- 1e-8
repeat{
  x1<- computeEstimate() #pvz kokia nors optimizavimo funkcija
  if(abs(x1-x0)<tol){
    break
  }else{
    x0<-x1
  }
}

#next comanda persoka prie sekancios interacijos
for(i in 1:100){
  if(i<=30){
    #skips the first 30
    next
  }
  print(i)
}

#sudeda skaicius kuriuos irasom i apacia add2(5,4) ir tt
add2<- function(x,y){
  x+y
}
add2(3,5)

#above10 tai skaiciai kurie yra >10. c- nurodo skaicius nuo 1 iki 20.
#above10(c) printina skaicius kurie atitinka salygas
above10<- function(x){
  use<- x>10
  x[use]
}
c<- seq(1:20)
above10(c)


#viskas panasiai kaip ir pries tai
above<- function(x,y=10){
  use<- x>y
  x[use]
}
above(c)
above(c,2)

#funkcija kuri apskaiciuoja stulpeliu vidurkius
column_mean<- function(y){
  nc<- ncol(y) #apskaiciuoja stulpeliu skaiciu
  means<- numeric(nc) 
  for(i in 1:nc){
    means[i]<- mean(y[,i])
  }
  means
}

library(datasets)
column_mean(airquality)

#kai kuriu stulpeliu neapskaiciavo nes buvo NA
#funkcija na.rm=TRUE ledizia apskaiciuoti pasalinus NA
column_mean<- function(y){
  nc<- ncol(y) 
  means<- numeric(nc) 
  for(i in 1:nc){
    means[i]<- mean(y[,i], na.rm=TRUE)
  }
  means
}
#praleidau pozicines vs leksines funkcijas bei datos funk - turbut nereik


#R 2.2

#loop funkcijos
#rasant scriptus for,while ir kiti loopai yra tinkami bet jeigu norima parasyti koda
#tiesiog konsoleje tada susiduriama su daug problemu

#lapply: loopina per list ir paleidzia funkcija kiekvienam elementui
#sapply: kaip ir lapply tik supaprastina rezultatus
#apply: taiko funkcija masyvo stulpeliams/eilutems
#tapply: taiko funkcija vektoriaus dalims
#mapply: multivariatine lapply funkcija

#lapply priima 3 argumentus (1) list objekta, (2) funkcija arba jos pavadinima
# (3) galimus funkcijos papildomus argumentus
#lapply visad grazina list klases objekta
x<- list(a=1:10, b=rnrom(10), c=seq(from=100, to=200, by=2))
lapply(x, mean)
#run if gerneuoja tolygiojo skirstinio atsitiktinius skaicius
lapply(x, runif, min=5, max=10)
# extra
z <- list(a=matrix(1:9, nrow=3, ncol = 3), 
          b=matrix(1:4, nrow = 2, ncol=2))
#elt yra anonimine funkcija. vietoj "elt" galetu buti "x", 
#arba "abrakadabra"
lapply(z, function(elt) elt[,1, drop=FALSE]) 

## sapply
#sapply` bando supaprastinti `lapply` rezultatus (jeigu imanoma)
# jeigu `lapply` grazintu list, kuriu kiekvienas elementas yra 1 ilgumo, tada `sapply` grazina vektoriu
# jeigu `lapply` grazintu list, kuriu kiekvienas elementas yra >1 ir vienodo ilgumo, tada `sapply` grazina matrica
# jeigu netinka pirma du variantai, grazina list kaip ir `lapply`
## sapply
set.seed(101)
x <-  list(a=1:10, 
           b=rnorm(50), 
           c=seq(from=100, to=200, by=2))
lapply(x, mean)
sapply(x, mean)


#apply funkcija taikoma matricu eilutems ir stulpeliams
#apply supaprastina loop naudojima
x <- matrix(1:4,2,2)
x
apply(x, 1, mean) # 1 - eiltuems
apply(x, 2, mean) # 2 - stulpeliams
apply(x, 1, sum) # 1 - eiltuems
apply(x, 2, sum) # 2 - stulpeliams

#jeigu norima apskaiciuoti dataframe/matricu eiluciu ar stulpeliu stumas/vidurkius
#galime naudoti jau supaprastintas funkcijas
# `rowSums=apply(x,1,sum)`
# `rowMeans=apply(x,1,mean)`
# `colSums=apply(x,2,sum)`
# `colMeans=apply(x,2,mean)

#norint pritaikyti apply funkcija daugiau dimensiju turinciam duomenu masyvui butina nurodyti vektoriu
# kurios dimensijos islaikomos
x<- array(data=rnorm(40), dim= c(2,2,10))
apply(x, c(1,2), mean)

#mapply taiko paraleliai(vienu metu) funkcija skirtingiems argumentams
list(rep(1,4), rep(2,3), rep(3,2), rep(4,1))
# tai galime supaprastinti su mapply
mapply(rep, 1:4, 4:1)
#funckija noise generuoja n atsitiktiniu normaliojo skirstinio skaiciu su vidurkiu mean ir standartiniu nuokrypiu sd
noise<- funtion(n,mean, sd){
  rnorm(n, mean, sd)
}
noise(4,1,2)

#mapply funkcija galima naudotis tam kad funkcija priimtu argumentus is vektoriu
noise<- funtion(n,mean, sd){
  rnorm(n, mean, sd)
}
mapply(noise, 1:4, 1:4, 0.1)

#tapply
#apply a function to each cell of a ragged array, that is to each (non empty) group of values given by a unique combinations of the levels of certain factors
#krc nesuprantu

#split padalina vektoriu arba kita objekta i grupes priklausomai nuo faktoriu arba faktoriu list

#distribucijos praleidau
#skirstiniai-
  #binominis atsitiktinis dydis(ad) modeliuoja viena kurio nors ivykio kurio sekmes tikmybe yra,p pasirodymu skaiciu atlikus n bandymu
  #posono - ad. modeliuojame kurio nors ivykio pasirodyma per fiksuota laiko itervala skaiciu(pvz:gautu per dienu sms zinuciu,auto avariju per men.ir tt)
  #tolygusis- tolygus intervalas nuo a iki b uzdaras ir norint aps jo vidurki paimi tuos 2 intervalo galus ir padalini is 2
  #normalusis- grafiko skirstinys yra varpo formos
  #eksponentinis ad. - daznai naudojamas laukimo laikams modeliujoti(pirmo kliento autobuso ir tt. kazko laukt vienu zodziu)

#3.1
#nustatyti adresa galima ir:
#setwd(../data) #taskiukai nustato per kiek folderiu perkelti adresa nuo dabartinio
#file.exists("file_name") testuoja ar egzistuoja failas arba direktorija
#dir.create("folder_name) - sukurioa direktorija 
#jeigu skriptas importuoja duomenis o duomenys turi buti patalpinti atskiroj direktorijoj:
if(!file.exists("data")){
  dir.create("data")
}
#paste() funkcija sujungia stringus/vectorius i viena character string
times<-"Times"; x<-3
paste("hello", "world",x,times,sep = "_")
#Duomenys is interneto
download.file(url,destfile,method,quiet=FALSE,mode="w",cacheOK = TRUE,extra=getOption("download.file.extra"))
#sitas tinka duomenu failu parsisiuntimui(.txt .csv, etc)
URL <- "http://atvira.sodra.lt/csv/lt-eur/apdraustieji_3_1.csv"
Sys.time()
Downloaddate <- format(Sys.time(), format = "%Y_%m_%d")
filename <- paste("./data/apdraustieji_3_1_",
                 Downloaddate,
                 ".csv",
                 sep="")
filename

download.file(url = URL,
              destfile = filename,
              method = "auto")


#dabar siunciames visa zip'a
URL <- "http://atvira.sodra.lt/downloads/lt-eur/apdraustuju_pajamu_analize.zip"
DownloadDate <- format(Sys.time(), format="%Y_%m_%d")
filename <- paste("./data/apdraustuju_pajamu_analize_",
                  DownloadDate,
                  ".zip", 
                  sep = "")
filename
download.file(
  url = URL, 
  destfile = filename,
  method = "auto")
unzip(filename,
      exdir = "data",
      files= "apdraustuju_pajamu_analize.csv")
GYV_PAJ<- read.csv("./data/apdraustuju_pajamu_analize.csv",
                   header=TRUE,
                   sep = ";",
                   stringsAsFactors = FALSE)
median(GYV_PAJ$Mënesio.pajamos)
list(lower_bound = 0.5*median(GYV_PAJ$Mënesio.pajamos),
     median = median(GYV_PAJ$Mënesio.pajamos),
     upper_bound = 2*median(GYV_PAJ$Mënesio.pajamos))
maximum<- 4000
hist(GYV_PAJ$M?nesio.pajamos[GYV_PAJ$M?nesio.pajamos<=maximum],
     main="Histogram of declared labor income",
     xlab="income grouped by 100",
     breaks=seq(0,maximum,100),
     xaxt="n")
axis(side=1, at=seq(0,maximum,100), labels = seq(0,maximum,100))

library("readdr")
install.packages("readr")

GYV_PAJ <- read.csv("./data/apdraustuju_pajamu_analize.csv",
                    header=TRUE,
                    sep=";",
                    stringsAsFactors = FALSE) 

median(GYV_PAJ$Mënesio.pajamos)
list(lower_bound = 0.5*median(GYV_PAJ$Mënesio.pajamos))
median= median(GYV_PAJ$Mënesio.pajamos)
upper_bound = 2*median(GYV_PAJ$Mënesio.pajamos)

maximum <- 4000
hist(GYV_PAJ$Mënesio.pajamos[GYV_PAJ$Mënesio.pajamos<=maximum],
     main="Histogram of declared labor income",
     xlab="Income grouped by 100",
     breaks= seq(0,maximum,100),
     xaxt="n")
axis(side= 1, at=seq(0,maximum,100), labels= seq(0,maximum,100))

df<- read.table("./data/apdraustuju_pajamu_analize.csv",
                header= TRUE,
                sep=";",
                fileEncoding= "ISO-8859-13",
                stringsAsFactors= FALSE)

#rsmdx sutvarko LSD failus
install.packages("rsdmx")
library(rsdmx)
#bbz nesamone

install.packages("eurostat")
library(eurostat)
nama_10_gdp <- get_eurostat("nama_10_gdp", stringsAsFactors = FALSE)

#duomenu inspektavimas
head(df,2)
tail(df,2)
str(df) #parodo struktura dataframe
summary(df) #issamesne struktura
quantile(df$Menesio.pajamos) #kazkoks kvantilius
quantile(df$Menesio.pajamos, probs = seq(0,1,0.1))
table(df$Lytis) #parodo kiek yra skirtingu reiksmiu tam stulepli
sum(is.na(df$Menesio.pajamos)) #parodo kiek yra NA reiksmiu
table(df$Lytis, df$Amžius) #susieja lenteles
 
#subsettingas

df<- mtcars
df[c(1:3),] #atspausdina pirmas 3 eilutes
df[c(1:3),c("mpg","cyl","hp")] #atspausdina pirmas tris eilutes su mph cyl ir hp reiksmem
df[df$mpg>=10& df$mpg<=20& df$hp>=250 | df$qsec<=16,] # | loginis OR

#which() grazina skaitini vektoriu priimdamas logini vektoriu
which(df$mpg<=14)
df[which(df$mpg<=14),1:3]
df[df$mpg<=14, 1:3]
head(df[order(df$mpg,df$hp),],3)
#sort() sortiruoja verktorius arba df su vienu kintamu
a <- c(1,NA,20,NA,40)
sort(a)
#order() su 2+ kintamais
head(df[order(df$mpg,df$hp),],3)

#nauju stulepliu kurimas
#ifelse() funkcija, jei taip tai taip, jei ne tai taip
mean(df$mpg)
df$loginis <- ifelse(df$mpg>=mean(df$mpg),"daugiau","mažiau")
head(df,3)

#cross tabs TAS GERAS PASIZIUREK DAR SUTVARKO FRAME
df <- read.csv("./data/apdraustieji_3_1_2019_06_17.csv",
               header=TRUE,
               sep=";",
               stringsAsFactors = FALSE)
df2 <- df[(df$Metai==2018 &
             df$Menesio.pajamos=="401-450 €"&
             df$Amžius=="Visos amžiaus grupes"&
             df$Lytis!="Visi apdraustieji"), ]


table <- xtabs(data=df2,
               Apdraustuju.skaicius ~ Lytis + Menuo,
               drop.unused.levels = TRUE)
table

barplot(table)

men <- c(Sausis=1,Vasaris=2,Kovas=3,
         Balandis=4,Geguže=5,Birželis=6,
         Liepa=7,Rugpjutis=8,Rugsejis=9,
         Spalis=10,Lapkritis=11,Gruodis=12)

df2$men <- men[df2$Menuo]
table <- xtabs(data=df2,
               Apdraustuju.skaicius ~ Lytis + men,
               drop.unused.levels = TRUE)
table
barplot(table)


#nuo cross tabu iki dplyr praleidau!!!!!!!!!!!!!!!
if(!require("dplyr")){
  install.packages("dplyr")
  library(dplyr)
}

library(eurostat)
nama_10_gdp <- get_eurostat("nama_10_gdp", stringsAsFactors = F, header= TRUE)

head(nama_10_gdp)
head(select(nama_10_gdp, 1:3)) #select leidzia pasirinkti variable's ar juos renamint ir tt
head(select(nama_10_gdp, unit, geo, values))
select(nama_10_gdp, -(1:3)) #isskyrus nuo pirmo iki 3

filter(nama_10_gdp, geo=="LT")

df <- filter(nama_10_gdp, geo=="LT" &
               na_item=="B1G" &
               unit=="CLV10_MEUR"&
               time>="2000-01-01")
plot(df$values, type="l")

df <- arrange(df, time) #apkeite data. buvo nuo virsaus i apacia db nuo apacios i virsu
plot(df$values, type="l")


df <- rename(df, rodiklis=na_item) # pakeisti pavadinima
head(df)

df <- mutate(df, mean=mean(df$values))
head(df)
plot(df$values, type="l")
lines(df$mean)

df <- filter(nama_10_gdp,
             na_item=="B1G" &
               unit=="CLV10_MEUR"&
               time>="2000-01-01")
df_g <- group_by(df, geo) #sugrupuoja pagal kazkoki kriteriju
summarise(df_g, mean=mean(values), median=median(values))
df
df_g
df <- nama_10_gdp %>%
  filter(geo=="LT" &
           na_item=="B1G" &
           unit=="CLV10_MEUR"&
           time>="2000-01-01")%>%
  select(time, values) %>%
  mutate(mean_LT=mean(values))%>% #sukuria arba transformina variables
  arrange(time) #apkeicia datos eiliskuma
plot(df$values, col="red", type="l")
lines(df$mean_LT, col="blue")

#----------------------------------------------------------
#4.1

# kazkokiareiksme %in% df$masinos parodo ar egzsituoja tas dalykas tame df
# | veikia kaip "arba"


# pasiziurek formules with, subset, par


par(mfrow=c(1,2))
with(airquality, plot(Ozone, Temp))
with(airquality, plot(Ozone, Wind ))

summary(ChickWeight)
x<- ChickWeight
with(subset(ChickWeight, Time==21), boxplot(weight~Diet))

with(subset(ChickWeight, Time==21), hist(weight))
with(subset(ChickWeight, Time==21), rug(weight)) #rug plottina pavienius elementus kaip bruksnelius ir leidzia suprasti ar deramai pasirinkti intervalai

######Base plotting funkcijos
#plot - sukuria pagrindini grafika
#lines - prideda linijas (vektorius)
#points - prideda taškus
#text - prideda teksta
#title - prideda anotacijas
#axis - prideda ašiu žymejimus, teksta

#####Base graphics parametrai
#par - gloablus parametrai
#bg - the background color
#mar - the margin size
#oma - the outer margin size
#mfrow - number of plots per row, column (filled row-wise)
#mfcol - number of plots per row, column (filled col-wise)



#plot(…):
#pch - the plotting symbol
#lty - the line type
#lwd - the line width
#col - color
#xlab - charackter string x-axis label
#ylab - charackter string y-axis label
#main - charackter string main label

###############################
## GGPLOT

library(ggplot2)
ggplot(mpg)+
  geom_point(aes(displ, hwy, color=drv)) # color= faktorius, drv duomenys suskaidomi pagal faktoriu ir nuspalvinami skirtingom spalvom

ggplot(mpg)+
  geom_point(aes(displ, hwy), color="red")# color="red" perdavus ne aes viduje visi taskai nudazomi raudonai

ggplot(mpg)+
  geom_point(aes(displ, hwy, color="red"))# perdavus viduje visi taskai nudazomi ir r galvoja kad tai kazkoks grupavimas

ggplot(mpg, aes(displ, hwy, color=drv))+ # jei aes() nekinta galime iskelti i ggplot dali
  
  geom_point()

# density plotai
if(!require("gridExtra")){
  install.packages("gridExtra")
  library(gridExtra)
  } #padeda sudeti kelis i viena
plot1 <- ggplot(mpg, aes(hwy))+geom_density()
plot2 <- ggplot(mpg, aes(hwy))+geom_histogram(binwidth = 2)
grid.arrange(plot1, plot2, ncol=2)

df<- data.frame(x=1:100, y=rnorm(100))
ggplot(df, aes(x,y))+theme_bw()+
  geom_line()+
  labs(x="Laikotapris", y="Ivertis", title= "Vidutine graža",
       subtitle = "Šaltinis Eurostat (nama_10_q)")

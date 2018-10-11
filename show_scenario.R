#########################################################################
#######Effect of antibiotic duration on hospitalised patients############
#########################################################################

rm(list=ls()) # Clean working environment

source("default_params.R")
source("model_binary.R")

####################6. Visualisation #####################

#present in one space
par(mfrow=c(4,2))

#plots
#antibiotics pixel
abx.mat.short<- abx.short[[3]]
abx.mat.long<- abx.long[[3]]
cols.a<-c('0'='grey', '1'='orange', '2'= 'orange3','3'= 'orange3')
abx_img.short<-image(1:nrow(abx.mat.short),1:ncol(abx.mat.short),abx.mat.short,col=cols.a, xlab='Time', ylab='Bed No.', main="Short duration of antibiotics")
abx_img.long<-image(1:nrow(abx.mat.long),1:ncol(abx.mat.long),abx.mat.long,col=cols.a, xlab='Time', ylab='Bed No.', main="Long duration of antibiotics")

#carriage pixel
orig_short <- colo_table_filled_short
dim.saved<-dim(orig_short)
colo_table_filled_short[colo_table_filled_short=="S" | colo_table_filled_short=="s"]<-1
colo_table_filled_short[colo_table_filled_short=="sr"|colo_table_filled_short=="Sr"|colo_table_filled_short=="sR"]<-2
colo_table_filled_short<-as.numeric(colo_table_filled_short)
dim(colo_table_filled_short)<-dim.saved

orig_long <- colo_table_filled_long
colo_table_filled_long[colo_table_filled_long=="S" | colo_table_filled_long=="s"]<-1
colo_table_filled_long[colo_table_filled_long=="sr"|colo_table_filled_long=="Sr"|colo_table_filled_long=="sR"]<-2
colo_table_filled_long<-as.numeric(colo_table_filled_long)
dim(colo_table_filled_long)<-dim.saved

cols<-c('chartreuse3', 'red')
col_img.short<-image(1:nrow(colo_table_filled_short),1:ncol(colo_table_filled_short), colo_table_filled_short, col=cols, xlab='Time', ylab='Bed No.')
col_img.long<-image(1:nrow(colo_table_filled_long),1:ncol(colo_table_filled_long),colo_table_filled_long,col=cols, xlab='Time', ylab='Bed No.')

#increase overtime
df.short<-as.data.frame(orig_short)
df.short$totalR<-rep(0, nrow(df.short))
for (i in 1:nrow(df.short)) {
    for (j in 1:ncol(df.short)) {
        if(df.short[i, j]=="sR"|df.short[i, j]=="sr"|df.short[i, j]=="Sr") {
            df.short$totalR[i] <- df.short$totalR[i] + 1
        }
    }
}

df.long<-as.data.frame(orig_long)
df.long$totalR<-rep(0, nrow(df.long))
for (i in 1:nrow(df.long)) {
    for (j in 1:ncol(df.long)) {
        if (df.long[i,j]=="sR"|df.long[i,j]=="sr"|df.long[i,j]=="Sr") {
            df.long$totalR[i] <- df.long$totalR[i] + 1
        }
    }
}

x <- c(1:n.days)
y.short <- df.short$totalR/n.bed
y.long <- df.long$totalR/n.bed
r.plot.short<-plot(x,y.short, type="l", ylim=c(0,1), xlab="Time", ylab="% patients carrying MDRO")
lines(x,y.long, col=2)
r.plot.long<-plot(x, y.long, type="l", ylim=c(0,1), xlab="Time", ylab="% patients carrying MDRO")

#total acquisitions 
df.short$newR <- rep(0, nrow(df.short))
df.long$newR <- rep(0, nrow(df.long))

for(i in 2:nrow(df.short)){
    for(j in 1:ncol(df.short)){
        if((df.short[i, j] == "sR"|df.short[i, j] == "sr"|df.short[i, j] == "Sr") & (df.short[i-1, j] == "S" | df.short[i-1, j] == "s")){
            df.short$newR[i] <- df.short$newR[i] + 1 
        }
    }
}
y.short<- df.short$newR

for(i in 2:nrow(df.long)){
    for(j in 1:ncol(df.long)){
        if((df.long[i, j] == "sR"|df.long[i, j] == "sr"|df.long[i, j] == "Sr") & (df.long[i-1, j] == "S" | df.long[i-1, j] == "s")){
            df.long$newR[i] <- df.long$newR[i] + 1 
        }
    }
}
y.long<- df.long$newR

acq.short<-plot(x,y.short, type='l', xlab="Time", ylab="Number of patients acquiring MDRO")
lines(x,y.long, col=2)
acq.long<-plot(x,y.long, type='l', xlab="Time", ylab="Number of patients acquiring MDRO")



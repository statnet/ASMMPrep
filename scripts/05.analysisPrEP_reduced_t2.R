
## CAMP Adol analysis file

library(mardham2)
library(EpiModelHPC)


s1 <- merge_simfiles(1, indir = "scenarios/adol/data/1to 99/")

s157 <- merge_simfiles(157, indir = "scenarios/adol/data/100 plus/")
s151 <- merge_simfiles(151, indir = "scenarios/adol/data/100 plus/")
s154 <- merge_simfiles(154, indir = "scenarios/adol/data/100 plus/")
s148 <- merge_simfiles(148, indir = "scenarios/adol/data/100 plus/")
s147 <- merge_simfiles(147, indir = "scenarios/adol/data/100 plus/")

s211 <- merge_simfiles(211, indir = "scenarios/adol/data/200 plus/")
s212 <- merge_simfiles(212, indir = "scenarios/adol/data/200 plus/")
s213 <- merge_simfiles(213, indir = "scenarios/adol/data/200 plus/")
s214 <- merge_simfiles(214, indir = "scenarios/adol/data/200 plus/")


#Base incidence.
incid.total.base<-rep(NA,100)

for (i in 1:100){
  
  incid.total.base[i]<-sum(s1$epi$incid[1041:2080,i])

}

incid.total.base<-mean(incid.total.base)




data<-c("s1","s157","s151","s154","s148","s147","s211", "s212", "s213",
        "s214") 



coverage<-rep(NA,length(data))
adherence<-rep(NA,length(data))
incid.total<-rep(NA,length(data))
incid.rate<-rep(NA,length(data))
incid.total.up95<-rep(NA,length(data))
incid.total.low95<-rep(NA,length(data))
incid.total.ly<-rep(NA,length(data))
incid.total.ly.up95<-rep(NA,length(data))
incid.total.ly.low95<-rep(NA,length(data))

prep.pt<-rep(NA,length(data))
prep.pt.up95<-rep(NA,length(data))
prep.pt.low95<-rep(NA,length(data))

person.time<-rep(NA,length(data))
person.time.deb<-rep(NA,length(data))

person.time.ly<-rep(NA,length(data))
person.time.ly.deb<-rep(NA,length(data))

prev.pop<-rep(NA,length(data))
prev.pop.up95<-rep(NA,length(data))
prev.pop.low95<-rep(NA,length(data))
prev.age18<-rep(NA,length(data))
prev.age18.up95<-rep(NA,length(data))
prev.age18.low95<-rep(NA,length(data))

PIA<-rep(NA,length(data))
PIA.low95<-rep(NA,length(data))
PIA.up95<-rep(NA,length(data))

NIA<-rep(NA,length(data))
NIA.low95<-rep(NA,length(data))
NIA.up95<-rep(NA,length(data))
NIA.deb<-rep(NA,length(data))
NIA.deb.low95<-rep(NA,length(data))
NIA.deb.up95<-rep(NA,length(data))

NNT<-rep(NA,length(data))
NNT.low95<-rep(NA,length(data))
NNT.up95<-rep(NA,length(data))

for (i in 1:length(data)){

  incid.total.temp<-rep(NA,100)
  incid.total.up95.temp<-rep(NA,100)
  incid.total.low95.temp<-rep(NA,100)
  incid.total.ly.temp<-rep(NA,100)
  incid.total.ly.up95.temp<-rep(NA,100)
  incid.total.ly.low95.temp<-rep(NA,100)
  prep.pt.temp<-rep(NA,100)
  prep.pt.up95.temp<-rep(NA,100)
  prep.pt.low95.temp<-rep(NA,100)
  person.time.temp<-rep(NA,100)
  person.time.deb.temp<-rep(NA,100)
  person.time.ly.temp<-rep(NA,100)
  person.time.ly.deb.temp<-rep(NA,100)
  coverage<-rep(NA,100)
  adherence<-rep(NA,100)
  
  prev.pop.temp<-rep(NA,100)
  prev.age18.temp<-rep(NA,100)
  
  PIA.temp<-rep(NA,100)
  NIA.temp<-rep(NA,100)
  NIA.deb.temp<-rep(NA,100)
  NNT.temp<-rep(NA,100)

  
  
    for (j in 1:100){  
      incid.total.temp[j]<-sum(get(data[i])$epi$incid[1041:2080,j])
      incid.total.ly.temp[j]<-sum(get(data[i])$epi$incid[2028:2080,j])
      prep.pt.temp[j]<-sum(get(data[i])$epi$prepCurr[1041:2080,j])
      person.time.temp[j]<-sum(get(data[i])$epi$num[1041:2080,j])
      person.time.deb.temp[j]<-sum(get(data[i])$epi$debuted[1041:2080,j])
      person.time.ly.temp[j]<-sum(get(data[i])$epi$num[2028:2080,j])
      person.time.ly.deb.temp[j]<-sum(get(data[i])$epi$debuted[2028:2080,j])
      prev.pop.temp[j]<-get(data[i])$epi$i.prev[2080,j]
      prev.age18.temp[j]<-get(data[i])$epi$i.prev.age6[2080,j]
      
      PIA.temp[j]<-(incid.total.base-incid.total.temp[j])/incid.total.base
      NIA.temp[j]<-((incid.total.base-incid.total.temp[j])/person.time.temp[j])*52*100000
      NIA.deb.temp[j]<-((incid.total.base-incid.total.temp[j])/person.time.deb.temp[j])*52*100000
      NNT.temp[j]<-(prep.pt.temp[j]/52)/(incid.total.base-incid.total.temp[j])

      

    }
  
  coverage<-get(data[i])$param$prep.coverage
  adherence<-get(data[i])$param$prep.class.prob
  incid.total[i]<-mean(incid.total.temp)
  incid.rate[i]<-mean(incid.total.temp)/1040
  incid.total.temp.ordered<-sort(incid.total.temp)
  incid.total.up95[i]<-incid.total.temp.ordered[length(incid.total.temp.ordered)*.97]
  incid.total.low95[i]<-incid.total.temp.ordered[length(incid.total.temp.ordered)*.03]
  
  incid.total.ly[i]<-mean(incid.total.ly.temp)
  incid.total.ly.temp.ordered<-sort(incid.total.ly.temp)
  incid.total.ly.up95[i]<-incid.total.ly.temp.ordered[length(incid.total.ly.temp.ordered)*.97]
  incid.total.ly.low95[i]<-incid.total.ly.temp.ordered[length(incid.total.ly.temp.ordered)*.03]
  
  prep.pt[i]<-mean(prep.pt.temp)
  prep.pt.temp.ordered<-sort(prep.pt.temp)
  prep.pt.up95[i]<-prep.pt.temp.ordered[length(prep.pt.temp.ordered)*.97]
  prep.pt.low95[i]<-prep.pt.temp.ordered[length(prep.pt.temp.ordered)*.03]
  
  prev.pop[i]<-mean(prev.pop.temp)
  prev.pop.temp.ordered<-sort(prev.pop.temp)
  prev.pop.up95[i]<-prev.pop.temp.ordered[length(prev.pop.temp.ordered)*.97]
  prev.pop.low95[i]<-prev.pop.temp.ordered[length(prev.pop.temp.ordered)*.03]
  
  prev.age18[i]<-mean(prev.age18.temp)
  prev.age18.temp.ordered<-sort(prev.age18.temp)
  prev.age18.up95[i]<-prev.age18.temp.ordered[length(prev.age18.temp.ordered)*.97]
  prev.age18.low95[i]<-prev.age18.temp.ordered[length(prev.age18.temp.ordered)*.03]
  
  
  person.time[i]<-mean(person.time.temp)
  person.time.deb[i]<-mean(person.time.deb.temp)
  
  person.time.ly[i]<-mean(person.time.ly.temp)
  person.time.ly.deb[i]<-mean(person.time.ly.deb.temp)
  
  PIA[i]<-mean(PIA.temp)
  PIA.temp<-sort(PIA.temp)
  PIA.low95[i]<-PIA.temp[length(PIA.temp)*.03]
  PIA.up95[i]<-PIA.temp[length(PIA.temp)*.97]
  
  NIA[i]<-mean(NIA.temp)
  NIA.temp<-sort(NIA.temp)
  NIA.low95[i]<-NIA.temp[length(NIA.temp)*.03]
  NIA.up95[i]<-NIA.temp[length(NIA.temp)*.97]
  
  NIA.deb[i]<-mean(NIA.deb.temp)
  NIA.deb.temp<-sort(NIA.deb.temp)
  NIA.deb.low95[i]<-NIA.deb.temp[length(NIA.deb.temp)*.03]
  NIA.deb.up95[i]<-NIA.deb.temp[length(NIA.deb.temp)*.97]
  
  NNT[i]<-mean(NNT.temp)
  NNT.temp<-sort(NNT.temp)
  NNT.low95[i]<-NNT.temp[length(NNT.temp)*.03]
  NNT.up95[i]<-NNT.temp[length(NNT.temp)*.97]
  

}

out.table2.supp<-cbind(PIA, PIA.low95, PIA.up95, NIA, NIA.low95, NIA.up95, NIA.deb, NIA.deb.low95, 
                       NIA.deb.up95, NNT, NNT.low95, NNT.up95, data)
out.table2.supp

library(xlsx) #load the package
write.xlsx(x = out.table2.supp, file = "scenarios/adol/out/table2.supp.xlsx",
           sheetName = "PrEP sims", row.names = FALSE)


##Make the output dataframe.

prep.out.matrix<-cbind(prev.age18,prev.age18.low95, prev.age18.up95, prev.pop, prev.pop.low95,prev.pop.up95,
                       incid.total,incid.total.low95,incid.total.up95,incid.total.ly,incid.total.ly.low95,incid.total.ly.up95,
                       prep.pt,prep.pt.low95,prep.pt.up95,person.time,person.time.deb,person.time.ly,person.time.ly.deb)

prep.out<-as.data.frame(prep.out.matrix)

#calulate incidence per 100 year at risk.
prep.out$incid.total.100pt.all<-NULL
prep.out$incid.total.low95.100pt.all<-NULL
prep.out$incid.total.up95.100pt.all<-NULL

prep.out$incid.total.100pt.deb<-NULL
prep.out$incid.total.low95.100pt.deb<-NULL
prep.out$incid.total.up95.100pt.deb<-NULL

prep.out$incid.total.ly.100pt.all<-NULL
prep.out$incid.total.ly.low95.100pt.all<-NULL
prep.out$incid.total.ly.up95.100pt.all<-NULL

prep.out$incid.total.ly.100pt.deb<-NULL
prep.out$incid.total.ly.low95.100pt.deb<-NULL
prep.out$incid.total.ly.up95.100pt.deb<-NULL


  prep.out$incid.total.100pt.all<-prep.out$incid.total/prep.out$person.time*52*100
  prep.out$incid.total.low95.100pt.all<-prep.out$incid.total.low95/prep.out$person.time*52*100
  prep.out$incid.total.up95.100pt.all<-prep.out$incid.total.up95/prep.out$person.time*52*100
  
  prep.out$incid.total.100pt.deb<-prep.out$incid.total/prep.out$person.time.deb*52*100
  prep.out$incid.total.low95.100pt.deb<-prep.out$incid.total.low95/prep.out$person.time.deb*52*100
  prep.out$incid.total.up95.100pt.deb<-prep.out$incid.total.up95/prep.out$person.time.deb*52*100
  
  prep.out$incid.total.ly.100pt.all<-prep.out$incid.total.ly/prep.out$person.time.ly*52*100
  prep.out$incid.total.ly.low95.100pt.all<-prep.out$incid.total.ly.low95/prep.out$person.time.ly*52*100
  prep.out$incid.total.ly.up95.100pt.all<-prep.out$incid.total.ly.up95/prep.out$person.time.ly*52*100
  
  prep.out$incid.total.ly.100pt.deb<-prep.out$incid.total.ly/prep.out$person.time.ly.deb*52*100
  prep.out$incid.total.ly.low95.100pt.deb<-prep.out$incid.total.ly.low95/prep.out$person.time.ly.deb*52*100
  prep.out$incid.total.ly.up95.100pt.deb<-prep.out$incid.total.ly.up95/prep.out$person.time.ly.deb*52*100
  
   base<-prep.out[1,]
  prep.out<-prep.out[2:length(data),]

  
#Number of infections averted per 100K person years at risk  (define at risk).
#Percent of infection averted.
#NNt persontime on prep / (1/NIA)
prep.out$NIA.all.pt<-rep(0,length(data)-1)
prep.out$NIA.deb.pt<-rep(0,length(data)-1)
prep.out$PIA<-rep(0,length(data)-1)
prep.out$NNT<-rep(0,length(data)-1)



#Number of infections averted per 100K person years at risk (in population and in the sexual marketplace).
for (i in 1:(length(prep.out$incid.total))){
    prep.out$NIA.all.pt[i]<-((incid.total.base-prep.out$incid.total[i])/prep.out$person.time[i])*52*100000
    prep.out$NIA.deb.pt[i]<-((incid.total.base-prep.out$incid.total[i])/prep.out$person.time.deb[i])*52*100000
    
    #NNT prep time over /gross number of infectionaverted..
    prep.out$NNT[i]<-(prep.out$prep.pt[i]/52)/(incid.total.base-prep.out$incid.total[i])
    }

#Percent of infection averted.
for (i in 1:(length(prep.out$incid.total))){
  prep.out$PIA[i]<-(incid.total.base-prep.out$incid.total[i])/incid.total.base
}



prep.out$incid.total
prep.out$sim<-data[2:length(data)]
prep.out

##Table 2.


table2.out<-cbind(prep.out$prev.age18, prep.out$prev.age18.low95, prep.out$prev.age18.up95, 
                  prep.out$prev.pop, prep.out$prev.pop.low95, prep.out$prev.pop.up95,
                  prep.out$incid.total.100pt.all, prep.out$incid.total.low95.100pt.all, prep.out$incid.total.up95.100pt.all,
                  prep.out$incid.total.100pt.deb, prep.out$incid.total.low95.100pt.deb, prep.out$incid.total.up95.100pt.deb,
                  prep.out$NIA.all.pt, prep.out$NIA.deb.pt, prep.out$PIA, prep.out$NNT,prep.out$prep.pt,prep.out$sim) 



library(xlsx) #load the package
write.xlsx(x = table2.out, file = "scenarios/adol/out/table2.xlsx",
           sheetName = "PrEP sims", row.names = FALSE)




save(prep.out, file = "scenarios/adol/out/sim.summary.output_t2.rda")



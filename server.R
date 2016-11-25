library(shiny)
shinyServer(function(input, output) {
      EMSaov.env$outputData<-NULL
      EMSaov.env$outANOVA<-NULL
      EMS.anova<-function(data.tot,Y.name,var.list,FixRan.list,nested.list=NULL,
                          model.level=NULL,n.table=NULL,approx.flag=FALSE,...){
        
        ## adjust the order of X variable for multi-level model
        
        if(!is.null(nested.list)){
          if(sum(!is.na(nested.list))!=0){
            nested.list<-unlist(lapply(nested.list,function(x){ 
              temp<-which(var.list==x);ifelse(length(temp)==0,NA,temp)}))
          }  
        }
        if(!is.null(model.level)){
          sort.id<-sort.list(model.level)
          nested.list<-var.list[nested.list[sort.id]]
          model.level<-model.level[sort.id]
          var.list<-var.list[sort.id]
          nested.list<-unlist(lapply(nested.list,function(x){ 
            temp<-which(var.list==x);ifelse(length(temp)==0,NA,temp)}))
          FixRan.list<-FixRan.list[sort.id]
          n.table[1:length(sort.id)]<-n.table[sort.id]
        }  
        if(is.null(nested.list)){
          nested.list<-rep(NA,length(var.list))
        }   
        if(is.null(n.table)){
          for(i in 1:length(var.list)){
            n.table<-c(n.table,length(table(data.tot[,var.list[i]])))
          }
          n.table<-c(n.table,mean(table(apply(data.tot[,var.list],1,
                                              function(x) paste(x,collapse="")))))
        }
        
        ## Change all X variables to factors
        
        data.tot<-data.tot[,c(var.list,Y.name)]
        for(i in var.list){
          data.tot[,i]<-factor(data.tot[,i])
        }
        
        ## design.M1
        
        n<-length(var.list)
        design.M1 <- NULL
        for(i in 1:n){
          design.M1<-rbind(design.M1,design.M1)
          temp1<-rep(c("",var.list[i]),each=2^(i-1))
          design.M1<-cbind(design.M1,temp1)
        }
        
        design.M1<-design.M1[-1,]
        
        ## Full model ANOVA
        
        model.F<-paste(Y.name,"~",paste(apply(design.M1,1,function(x) 
          paste(paste(x[x!=""],collapse="*"))),collapse="+"))
        model.id<-c(apply(design.M1,1,function(x) 
          paste(paste(x[x!=""],collapse=":"))),"Residuals")
        options(warn=-1)
        SS.table<-stats::anova(stats::lm(eval(model.F),
                                         data = data.tot))[model.id,1:2]
        options(warn=0)
        ## treat nested
        
        colnames(design.M1)<-var.list
        nest.id<-which(!is.na(nested.list))
        
        if(length(nest.id)>0){
          for(i in 1:length(nest.id)){
            temp.list<-var.list[apply(design.M1[,1:n],1,function(x) 
              ifelse(sum(x==var.list[nest.id[i]])==0,
                     NA,nested.list[nest.id[i]]))]
            del.list<-which(apply(design.M1[,1:n],1,function(x) 
              sum(x==var.list[nest.id[i]])*
                sum(x==var.list[nested.list[nest.id[i]]]))==1)
            for(k in 1:length(del.list)){
              comb.id<-del.list[k]
              temp.k<-design.M1[comb.id,]
              temp.k<-temp.k[temp.k!="" & temp.k!=var.list[nested.list[nest.id[i]]]]
              temp.k<-paste(temp.k,collapse="")
              comb.id<-c(comb.id,which(apply(design.M1,1,
                                             function(x) paste(x,sep="",collapse=""))==temp.k))
              SS.temp<-apply(SS.table[comb.id,],2,sum)
              SS.table[comb.id[length(comb.id)],]<-SS.temp
            }
            design.M1<-cbind(design.M1,temp.list)
            colnames(design.M1)<-c(colnames(design.M1)[-ncol(design.M1)],"nested")
            design.M1<-design.M1[-del.list,]
            SS.table<-SS.table[-del.list,]
            
            ## nested-nested-...
            
            flag<-TRUE
            id.t<-nest.id[i]
            while(flag){
              temp.c<-nested.list[nested.list[id.t]]
              if(is.na(temp.c)){
                flag<-FALSE
              }else{
                del.list<-which(apply(design.M1,1,function(x) 
                  sum(x[1:n]==var.list[temp.c])*!is.na(x[n+i]))==1)
                design.M1<-design.M1[-del.list,]    
                SS.temp<-apply(SS.table[del.list,],2,sum)
                design.M1[which(design.M1[,n+i] ==
                                  var.list[nested.list[nest.id[i]]]),
                          n+temp.c]<-var.list[temp.c]
                sel.id<-which(design.M1[,n+i] ==
                                var.list[nested.list[nest.id[i]]])
                SS.table[sel.id,] <-SS.table[sel.id,]+SS.temp
                SS.table<-SS.table[-del.list,]
              }
              id.t<-nested.list[id.t]
            }
          } 
        }  
        
        ## EMS.table
        
        design.M1[is.na(design.M1)]<-""
        out<-apply(design.M1,1,function(x) 
          ifelse(paste(x[-(1:n)],collapse="")!="",
                 paste(paste(x[1:n][x[1:n]!=""],collapse=":"),
                       "(",paste(x[-(1:n)][x[-(1:n)]!=""],collapse=","),
                       ")",sep=""),
                 paste(x[1:n][x[1:n]!=""],collapse=":")))
        rownames(SS.table)[-nrow(SS.table)]<-out
        EMS.table<-matrix(0,ncol=length(var.list)+1,nrow=length(out)+1)
        colnames(EMS.table)<-c(var.list,"Error")
        rownames(EMS.table)<-c(out,"Error")
        n.EMS<-nrow(EMS.table)
        p.EMS<-ncol(EMS.table)
        EMS.table[,p.EMS]<-n.table[p.EMS]
        EMS.table[n.EMS,]<-1
        temp<-design.M1[,1:length(var.list)]
        temp.nest<-design.M1[,-c(1:length(var.list)),drop=FALSE]
        temp[temp==""]<-NA
        for(i in 1:ncol(temp)){
          if(sum(temp[,i]==var.list[i],na.rm=TRUE)!=0){
            id.t<-which(is.na(temp[,i]))
            EMS.table[id.t,i]<-n.table[i]
            if(FixRan.list[i]=="R")  EMS.table[-c(id.t,n.EMS),i]<-1
          }else{
            sel.id<-which(!is.na(temp[,i]))
            EMS.table[sel.id,which(var.list==temp[sel.id,i][1])]<-1
          }   
          if(length(nest.id)>0){
            for(k in 1:ncol(temp.nest))
              EMS.table[which(temp.nest[,k]==var.list[i]),i]<-1     
          }
        }    
        
        ## EMS  
        
        temp.t<-design.M1[,1:length(var.list)]
        EMS<-NULL
        n.E<-nrow(EMS.table)
        id.keep<-NULL 
        hid.flag<-NULL
        for(i in n.E:1){
          if(i!=n.E){
            sel.id<-temp.t[i,,drop=FALSE]
            if(length(nest.id)>0){
              tt<-temp.nest[i,]
              id.keep<-NULL
              for(l in 1:length(tt))
                id.keep<-c(id.keep,which(names(hid.flag)==tt[l]))
            }          
            hid.flag<-rep(TRUE,ncol(EMS.table))
            names(hid.flag)<-colnames(EMS.table)
            for(j in 1:length(var.list)){
              hid.flag[which(names(hid.flag)==sel.id[j])]<-FALSE
            }
            pick.id<-design.M1[i,]
            pick.id<-pick.id[pick.id!=""]
            temp<-apply(design.M1,1,function(x) {
              keep.t<-TRUE; 
              for(i in 1:length(pick.id)) 
                keep.t<-keep.t*(sum(x==pick.id[i])!=0)
              return(keep.t)})
            hid.flag[id.keep]<-TRUE
            temp.T<-apply(EMS.table[,hid.flag,drop=FALSE],1,prod)
            temp.T.1<-temp.T[temp.T!=0]
            temp.T.1[length(temp.T.1)]<-""
            name.temp.T<-names(temp.T.1)
            temp<-c(temp,1)[temp.T!=0]
            nn<-length(temp.T.1)    
            temp.T.1<-temp.T.1[nn:1] 
            name.temp.T<-name.temp.T[nn:1] 
            temp<-temp[nn:1] 
            temp.EMS<-paste(temp.T.1[temp==1],name.temp.T[temp==1],
                            sep="",collapse="+")
            
          }else{
            temp<-c(rep(0,ncol(EMS.table)-1),1)
            temp.EMS<-"Error"
          }
          EMS<-cbind(temp.EMS,EMS)
        }
        
        ## model level
        
        if(!is.null(model.level)){
          level.list<-sort(unique(model.level))
          n.L<-length(level.list)
          Model.level<-rep(level.list[n.L],nrow(SS.table)-1)    
          temp.flag<-rep(TRUE,length(Model.level))
          for(i in n.L:1){
            i.id<-which(model.level==i)
            for(k in i.id){
              Model.level[which((design.M1[,k]!="")*temp.flag==1)]<-i
              temp.flag[design.M1[,k]!=""]<-FALSE
            }    
          }
          Model.level<-c(Model.level,max(Model.level))
        }else{
          Model.level<-NULL
        }   
        n.t<-nrow(SS.table)
        flag.zero.MSE<-FALSE  
        
        if(SS.table[n.t,2]==0){
          SS.table[n.t,1:2]<-SS.table[n.t-1,1:2] 
          temp.name<-rownames(SS.table)[n.t]
          SS.table<-SS.table[-n.t,]
          rownames(SS.table)[n.t-1]<-temp.name
          t.EMS<-lapply(EMS,function(x) strsplit(x,"[+]")[[1]])
          del.list<-t.EMS[[n.t-1]][-1]
          for(i in 1:length(t.EMS)){
            keep.id<-NULL
            for(j in 1:length(del.list))
              keep.id<-c(keep.id,which(t.EMS[[i]]==del.list[j]))
            if(length(keep.id)!=0)
              t.EMS[[i]]<-t.EMS[[i]][-keep.id]
          }
          EMS<-unlist(lapply(t.EMS,function(x) paste(x,sep="",collapse="+")))
          EMS[n.t-1]<-EMS[n.t]
          EMS<-EMS[1:(n.t-1)]
          flag.zero.MSE<-TRUE
          Model.level<-Model.level[1:(n.t-1)]
        }
        ## Caldulate MS, F (approx.F), P-value, sig,
        
        SS.table[,3]<-SS.table[,2]/SS.table[,1]
        split.EMS<-lapply(EMS,function(x) strsplit(x,"[+]")[[1]])
        F.value<-NULL
        P.value<-NULL
        Signif<-NULL
        for(i in 1:nrow(SS.table)){
          n.SE<-length(split.EMS[[i]])
          SS.temp<-paste(split.EMS[[i]][-n.SE],collapse="+")
          if(sum(EMS==SS.temp)!=0){
            F.temp<-SS.table[i,3]/SS.table[which(EMS==SS.temp),3]
            pValue.temp<- 1-stats::pf(F.temp,SS.table[i,1],
                                      SS.table[which(EMS==SS.temp),1])
          }else if(i!=nrow(SS.table) & approx.flag){
            test.EMS<-split.EMS[[i]]
            Appr.result<-Approx.F(SS.table=data.frame(SS.table,EMS=c(EMS)),approx.id=i)
            F.temp<-Appr.result$Appr.F
            pValue.temp<-Appr.result$Appr.Pvalue
          } else{
            F.temp<-NA
            pValue.temp<-NA
          }
          
          if(!is.na(pValue.temp)){
            if(pValue.temp<=0.001){
              Signif.temp <- "***"
            }else if(pValue.temp<=0.01){
              Signif.temp <- "**"
            }else if(pValue.temp<=0.05){
              Signif.temp <- "*"
            }else if(pValue.temp<=0.1){
              Signif.temp <- "."
            }else{
              Signif.temp <- ""
            }
            pValue.temp <- ifelse(round(pValue.temp,4)<0.0001,
                                  "<0.0001",round(pValue.temp,4))
            F.temp <- round(F.temp,4)
          }else{
            Signif.temp <- ""
            pValue.temp <- ""
            F.temp<-""
          }
          F.value<-c(F.value,F.temp)
          P.value<-c(P.value,pValue.temp)
          Signif<-c(Signif,Signif.temp)    
        }
        
        SS.table.t<-cbind(SS.table[,1],
                          round(SS.table[,2],4),
                          round(SS.table[,3],4))
        colnames(SS.table.t)<-c("Df","SS","MS")
        if(!is.null(Model.level)){
          tot.result<-data.frame(SS.table.t,Fvalue=F.value,Pvalue=P.value,
                                 Sig=Signif,Model.Level=Model.level,EMS=matrix(EMS))    
        }else{
          tot.result<-data.frame(SS.table.t,Fvalue=F.value,Pvalue=P.value,
                                 Sig=Signif,EMS=matrix(EMS))   
        }
        rownames(tot.result)<-rownames(SS.table)  
        return(tot.result)
      }
      
      Approx.F<-function(SS.table,approx.id,...){
        EMS<-as.character(SS.table$EMS)
        split.EMS<- lapply(EMS,function(x) strsplit(x,"[+]")[[1]])
        split.EMS.last<-lapply(split.EMS,function(x) return(x[length(x)]))  
        test.EMS<-split.EMS[[approx.id]]
        n.SE<-length(test.EMS)
        TEMP.EMS<-test.EMS[-n.SE]
        keep.id<-NULL
        keep.var<-NULL
        for(kk in 2:length(TEMP.EMS)){
          keep.id<-c(keep.id,which(split.EMS.last==TEMP.EMS[kk]))
          keep.var<-c(keep.var,TEMP.EMS[kk])  
        }
        keep.var<-keep.var[keep.id!=1]
        TEMP.EMS<-unlist(split.EMS[keep.id])
        TEMP.EMS<-TEMP.EMS[TEMP.EMS!="Error"]
        den.id<-names(table(TEMP.EMS))[table(TEMP.EMS)==1]
        ms.num<-SS.table[approx.id,3]
        ms.den<-0
        df.num<-SS.table[approx.id,3]^2/SS.table[approx.id,1]
        df.den<-0
        
        for(kk in 1:length(keep.var)){
          if(sum(keep.var[kk]==den.id)==1){
            id.i<-which(split.EMS.last==keep.var[kk])
            ms.den<-ms.den+SS.table[id.i,3]
            df.den<-df.den+SS.table[id.i,3]^2/SS.table[id.i,1]
          }else{  
            id.i<-which(split.EMS.last==keep.var[kk])
            ms.num<-ms.num+SS.table[id.i,3]
            df.num<-df.num+SS.table[id.i,3]^2/SS.table[id.i,1]    
          }
        }
        Appr.F<-ms.num/ms.den
        Appr.F.df1<-ms.num^2/df.num
        Appr.F.df2<-ms.den^2/df.den
        Appr.Pvalue<-1-stats::pf(Appr.F,Appr.F.df1,Appr.F.df2)
        return(list(Appr.F=Appr.F,Appr.Pvalue=Appr.Pvalue))
      }  
      
      
      PooledANOVA<-function(SS.table,del.ID,...){
        temp.SS<-SS.table[,1:2]
        temp.EMS<-as.character(SS.table$EMS)
        Model.level<-SS.table$Model.Level
        temp.ID<-del.ID[del.ID!="Residuals"]
        temp.ID<-unlist(lapply(temp.ID,function(x) which(rownames(temp.SS)==x)))
        temp.EMS<-as.character(temp.EMS)
        temp.SS[nrow(temp.SS),]<-apply(temp.SS[del.ID,],2,sum)
        temp.SS<-temp.SS[-temp.ID,]
        Model.level<-Model.level[-temp.ID]
        
        temp.SS[,3]<-temp.SS[,2]/temp.SS[,1]
        temp.split.EMS<-lapply(temp.EMS,function(x) {
          temp1<-strsplit(x,"[+]")[[1]]
          for(i in 1:length(temp.ID)){
            t.id<-grep(del.ID[i],temp1)
            if(length(t.id)!=0)
              temp1<-temp1[-t.id]
          }
          return(temp1)})
        
        temp.split.EMS<-temp.split.EMS[-temp.ID]    
        EMS.t<-lapply(temp.split.EMS,function(x) paste(x,sep="",collapse="+"))
        
        F.value<-NULL
        P.value<-NULL
        Signif<-NULL
        
        for(i in 1:nrow(temp.SS)){
          n.SE<-length(temp.split.EMS[[i]])
          SS.temp<-paste(temp.split.EMS[[i]][-n.SE],collapse="+")
          test.EMS<-temp.split.EMS[[i]]
          if(sum(temp.EMS==SS.temp)!=0){
            F.temp<-temp.SS[i,3]/temp.SS[which(EMS.t==SS.temp),3]
            pValue.temp<- 1-stats::pf(F.temp,temp.SS[i,1],
                                      temp.SS[which(EMS.t==SS.temp),1])
          } else if(i!=nrow(temp.SS)&length(test.EMS)!=1){
            Appr.result<-Approx.F(data.frame(temp.SS,EMS=unlist(EMS.t)),1)
            F.temp<-Appr.result$Appr.F
            pValue.temp<-Appr.result$Appr.Pvalue
          } else{
            F.temp<-NA
            pValue.temp<-NA
          }
          if(!is.na(pValue.temp)){
            if(pValue.temp<=0.001){
              Signif.temp <- "***"
            }else if(pValue.temp<=0.01){
              Signif.temp <- "**"
            }else if(pValue.temp<=0.05){
              Signif.temp <- "*"
            }else if(pValue.temp<=0.1){
              Signif.temp <- "."
            }else{
              Signif.temp <- ""
            }
            pValue.temp <- ifelse(round(pValue.temp,4)<0.0001,
                                  "<0.0001",round(pValue.temp,4))
            F.temp <- round(F.temp,4)
          }else{
            Signif.temp <- ""
            pValue.temp <- ""
            F.temp<-""
          }
          F.value<-c(F.value,F.temp)
          P.value<-c(P.value,pValue.temp)
          Signif<-c(Signif,Signif.temp)    
        }
        
        SS.table.t<-cbind(temp.SS[,1],
                          round(temp.SS[,2],4),
                          round(temp.SS[,3],4))
        colnames(SS.table.t)<-c("Df","SS","MS")
        EMS.t<-as.character(EMS.t)
        if(!is.null(Model.level)){
          tot.result<-data.frame(SS.table.t,Fvalue=F.value,Pvalue=P.value,
                                 Sig=Signif,Model.Level=Model.level,EMS=matrix(EMS.t))    
        }else{
          tot.result<-data.frame(SS.table.t,Fvalue=F.value,Pvalue=P.value,
                                 Sig=Signif,EMS=matrix(EMS.t))   
        }
        rownames(tot.result)<-rownames(temp.SS) 
        return(tot.result)
      }
      
      
      Dataset<-shiny::reactive({
        if(is.null(input$outputfile)){
          return(data.frame())
        }
       EMSaov.env$outputData<-data.frame(do.call("read.csv",
                                        list(input$outputfile$datapath))) 
        return(EMSaov.env$outputData)
      })
     
      output$choose_Yvar<-shiny::renderUI({
        if(is.null(input$outputfile))
          return()
        if(identical(Dataset(),'')||identical(Dataset(),data.frame())) 
          return(NULL)      
       EMSaov.env$outputData<-Dataset()
       EMSaov.env$NUM<-dim(EMSaov.env$outputData)[2] #Num of all variable  ##
       EMSaov.env$Class<-sapply(apply(EMSaov.env$outputData,2,unique),length)  ##  
        EMSaov.env$Colnames<-colnames(EMSaov.env$outputData)##
        shiny::selectInput("Yvar",label="Y variable",c("",EMSaov.env$Colnames))    
      })
     
      output$choose_Xvar<-shiny::renderUI({
        if(is.null(input$outputfile))
          return()    
        if(is.null(input$outputfile)|is.null(EMSaov.env$outputData)){
          choice.temp<-c(" "," ") 
        }else{
          choice.temp<-c(EMSaov.env$Colnames)
        }      
        shiny::checkboxGroupInput("Xvar","X variable",choices=choice.temp) 
      })  

      output$choose_type<-shiny::renderUI({
        if(is.null(input$outputfile))
          return()    
        if(is.null(input$outputfile)|is.null(EMSaov.env$outputData)){
          choice.temp<-c(" "," ")
        }else{
          choice.temp<-c(EMSaov.env$Colnames)
        } 
        shiny::checkboxGroupInput("type","Random Effect",choices=choice.temp) 
      }) 

      makenumericButton<-function(n){
        if(n==1){
          shiny::numericInput(paste0("level",n),
                      label=paste0("[# of categories] ",EMSaov.env$Colnames[n]),
                      value=EMSaov.env$Class[n])
        }else{
          shiny::numericInput(paste0("level",n),label=EMSaov.env$Colnames[n],
                              value=EMSaov.env$Class[n])
        }
      }
      WidgetVector<-shiny::reactive({lapply(X=1:EMSaov.env$NUM,
                                            FUN=makenumericButton)})
      output$choose_level<-shiny::renderUI({
        if(is.null(input$outputfile)|is.null(EMSaov.env$outputData)){
          return()  
        }else{
          shiny::tagList(WidgetVector())
        }
      }) 
     
      makeselectButton<-function(n){
        if(n==1){
          shiny::selectInput(paste0("nested",n),
                     label=paste0("[nested]\n ",EMSaov.env$Colnames[n]),
                     c("None",EMSaov.env$Colnames))
        }else{
          shiny::selectInput(paste0("nested",n),
                             label=EMSaov.env$Colnames[n],
                             c("None",EMSaov.env$Colnames))
        }
      }
     
      WidgetVector2<-shiny::reactive({lapply(X=1:EMSaov.env$NUM,
                                             FUN=makeselectButton)})
      output$choose_nested<-shiny::renderUI({
        if(is.null(input$outputfile)| is.null(EMSaov.env$outputData)){
          return()  
        }else{
          shiny::tagList(WidgetVector2())
        }
      }) 
     
      makenumericButton2<-function(n){
        if(n==1){
          shiny::numericInput(paste0("split",n),
                      label=paste0("[model level] ",EMSaov.env$Colnames[n]),
                      value=1)
        }else{
          shiny::numericInput(paste0("split",n),
                              label=EMSaov.env$Colnames[n],value=1)
        }
      }
     
      WidgetVector3<-shiny::reactive({lapply(X=1:EMSaov.env$NUM,
                                             FUN=makenumericButton2)})
      output$choose_split<-shiny::renderUI({
        if(is.null(input$outputfile)|is.null(EMSaov.env$outputData)){
          return()  
        }else{
          shiny::tagList(WidgetVector3())
        }
      }) 
     
      output$EDA1<-shiny::renderPlot({
        if(is.null(input$outputfile)|is.null(EMSaov.env$outputData)| 
           is.null(input$Xvar)|is.null(input$Yvar)){
          return()
        }else{
          X<-EMSaov.env$outputData[,input$Xvar]
          Y<-EMSaov.env$outputData[,input$Yvar]
          p<-length(input$Xvar)
          r<-ceiling(sqrt(p))
          graphics::par(mfrow=c(1,p))
          for(i in 1:p){
            graphics::plot(Y~factor(X[,i]),xlab=input$Xvar[i],ylab=input$Yvar)
            graphics::points(1:length(table(X[,i])),tapply(Y,X[,i],mean),
                             col=2,pch=16,cex=1.5)
          }  
        }
      })
     
      output$EDA2<-shiny::renderPlot({
        if(is.null(input$outputfile)|is.null(EMSaov.env$outputData)| 
           is.null(input$Xvar)|is.null(input$Yvar)){
          return()
        }else{
          X<-EMSaov.env$outputData[,input$Xvar]
          Y<-EMSaov.env$outputData[,input$Yvar]
          p<-length(input$Xvar)
          r<-ceiling(sqrt(p*(p-1)/2))
          graphics::par(mfrow=c(r,r))
          for(i in 1:(p-1)){
            for(j in (i+1):p){
              temp.group<-as.numeric(X[,j])
              r<-length(table(X[,i]))
              graphics::matplot(c(-0.5,r),range(Y),type="n",
                     xlab=input$Xvar[i],ylab=input$Yvar,
                     main=paste(input$Xvar[i],"*", input$Xvar[j]))
              temp.table<-names(table(X[,j]))
              for(k in 1:length(temp.table)){
                graphics::lines(1:length(table(X[temp.group==k,i])),
                  tapply(Y[temp.group==k],X[temp.group==k,i],mean),lty=k,col=k)
              }
              graphics::legend(-0.5,max(Y),temp.table,lty=1:r,col=1:r,
                               title=input$Xvar[j])
            }
          }
        }
      })
     
      output$result1<-shiny::renderTable({
        if(is.null(input$outputfile)| is.null(EMSaov.env$outputData) | 
           is.null(input$Xvar)| is.null(input$Yvar)){
          return()
        }else{
          X<-EMSaov.env$outputData[,input$Xvar]
          Y<-EMSaov.env$outputData[,input$Yvar]
          for(i in 1:EMSaov.env$NUM){
           EMSaov.env$Class[i]<-input[[paste0("level",i)]]
          }  #inputEMSaov.env$Class   
          level<-EMSaov.env$Class[c(input$Xvar)]
          level<-c(level,mean(table(X)))

          Type<-matrix("F",nrow=length(input$Xvar))
          rownames(Type)<-input$Xvar
          Type[input$type,]<-"R"
          type<-c(Type)
         
          nested<-NULL
          for(i in 1:EMSaov.env$NUM){
            nest<-input[[paste0("nested",i)]]
            if(is.null(nest)){
              nested[i]<-""
            }else{
              nested[i]<-nest
            }
          }
          names(nested)<-EMSaov.env$Colnames
          nested<-nested[input$Xvar]
          n<-length(input$Xvar)
         
         #split
          split<-NULL
          for(i in 1:EMSaov.env$NUM)
            split[i]<-input[[paste0("split",i)]]
          names(split)<-EMSaov.env$Colnames
         
          split<-split[c(input$Xvar)]
          split<-split[!is.na(split)]
          var.list<-input$Xvar      
          nest.temp<-rep(NA,length(nested))
          for(i in 1:length(nested))
            nest.temp[i]<-ifelse(nested[i]=="",NA,which(var.list==nested[i]))
          nested<-nest.temp
          if(sum(split==1)==length(split)) 
            split<-NULL      
          data.tot<-EMSaov.env$outputData[,c(input$Xvar,input$Yvar)]
          out<- EMS.anova(data.tot=data.tot,
                                 Y.name=input$Yvar,
                                 var.list=input$Xvar,
                                 FixRan.list=type,                        
                                 nested.list=input$Xvar[nested],
                                 model.level=split)
        }      
      })
     
      output$result2<-shiny::renderTable({
        if(is.null(input$outputfile)|is.null(EMSaov.env$outputData)| 
           is.null(input$Xvar)| is.null(input$Yvar)){
          return()
        }else{
          X<-EMSaov.env$outputData[,input$Xvar]
          Y<-EMSaov.env$outputData[,input$Yvar]
         
          for(i in 1:EMSaov.env$NUM){
           EMSaov.env$Class[i]<-input[[paste0("level",i)]]
          }  #inputEMSaov.env$Class     
         
          level<-EMSaov.env$Class[c(input$Xvar)]
          level<-c(level,mean(table(X)))
          Type<-matrix("F",nrow=length(input$Xvar))
          rownames(Type)<-input$Xvar
          Type[input$type,]<-"R"
          type<-c(Type)
         
         #nested   
          nested<-NULL
          for(i in 1:EMSaov.env$NUM){
            nest<-input[[paste0("nested",i)]]
            if(is.null(nest)){
              nested[i]<-""
            }else{
              nested[i]<-nest
            }
          }
          names(nested)<-EMSaov.env$Colnames
          nested<-nested[input$Xvar]
          n<-length(input$Xvar)
          split<-NULL
          for(i in 1:EMSaov.env$NUM)
            split[i]<-input[[paste0("split",i)]]
          names(split)<-EMSaov.env$Colnames
          split<-split[c(input$Xvar)]
          split<-split[!is.na(split)]
          var.list<-input$Xvar      
          nest.temp<-rep(NA,length(nested))
          for(i in 1:length(nested))
            nest.temp[i]<-ifelse(nested[i]=="None",
                                 NA,which(var.list==nested[i]))
          nested<-nest.temp
          if(sum(split==1)==length(split)) 
            split<-NULL      
          data.tot<-EMSaov.env$outputData[,c(input$Xvar,input$Yvar)]
          out<- EMS.anova(data.tot=data.tot,
                                 Y.name=input$Yvar,
                                 var.list=input$Xvar,
                                 FixRan.list=type,                        
                                 nested.list=input$Xvar[nested],
                                 model.level=split,
                                 approx.flag=TRUE)
          EMSaov.env$outANOVA<-out
        }      
      })
     
      output$choose_ANOVA <-  shiny::renderUI({
        if(is.null(input$outputfile))
          return()
        if(identical(Dataset(),'')||identical(Dataset(),data.frame())||
           is.null(EMSaov.env$outANOVA)) 
          return(NULL)  
        Rnames<-rownames(EMSaov.env$outANOVA)
        print(Rnames)
        shiny::checkboxGroupInput("ANOVA","Combine ANOVA table",choices=Rnames)  
      })
     
      output$result3<-shiny::renderTable({
        if(is.null(input$outputfile)|is.null(EMSaov.env$outputData)| 
           is.null(input$Xvar)|is.null(input$Yvar)|is.null(input$ANOVA)){
          return()
        }else{
          sel.id<-NULL
          temp.input<-unique(c(input$ANOVA,"Residuals"))
          for(i in temp.input)
            sel.id<-c(sel.id,which(rownames(EMSaov.env$outANOVA)==i))
          if(length(sel.id)>1){
            temp.SS<-EMSaov.env$outANOVA[,1:2]
            temp.SS$Df<-as.numeric(as.character(temp.SS$Df))
            temp.SS$SS<-as.numeric(as.character(temp.SS$SS))
            Residuals<-apply(temp.SS[sel.id,],2,sum)
            temp.SS<-rbind(temp.SS[-sel.id,],Residuals)
            rownames(temp.SS)[nrow(temp.SS)]<-"Residuals"
            temp.EMS<-c(as.character(EMSaov.env$outANOVA$EMS)[-sel.id],
               as.character(EMSaov.env$outANOVA$EMS)[nrow(EMSaov.env$outANOVA)])
            del.ID<-temp.input
            Model.level<-EMSaov.env$outANOVA$Model.Level
            Model.level<-c(Model.level[-sel.id],
                           Model.level[length(Model.level)])
            out<- PooledANOVA(EMSaov.env$outANOVA,del.ID)
          }else{
            out<-EMSaov.env$outANOVA 
          }
        }      
      })
})
      


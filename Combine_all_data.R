# Coursework_Team_12
## Insert the data and make into a table
install.packages("janitor")
library(janitor)
install.packages("dplyr")
library(dplyr)

getwd() #"/Users/miasmacbook/Desktop/Group12/RScript"
setwd(dir = "Desktop//Group12/data")
test_1 <- read.csv("data/0a00m3t298175av.csv", header =F)
test_1 <- t(test_1)
test_1 <- as.data.frame(test_1)
test_1 <- row_to_names(test_1, row_number = 1)
sort_t1 <- test_1%>%select(order(colnames(.)))

#get start with real data 
#list the file in the folder
files <- list.files(path="./data", pattern="*.csv", full.names = TRUE, recursive=FALSE)
#cuz transport the data will loss the files (i don't know why, but everytime i did it, my files will be 53 only)
#therefore leave the t() at the end
#then need to grab the first column out for every file, cuz it is the column name in the final doc
first <- as.data.frame(read.csv(files[1], header = F))
#sort the table by the first column (a-z)
first <- first[order(first[[1]]),]
#save the first column as the samecol, then can use it individually 
samecol <- first[,1]
samecol
#Here is another checking processing right here
#checking can i add the samecol as the first column before transport
first_fin <- first[,-1]
first_fin
test <- data.frame(col=samecol,stringsAsFactors = FALSE)
test <- cbind(test, first_fin)
test <- t(test)
test <- row_to_names(test, row_number = 1)

#get result as i want
#start the processes for whole folder 
#set up a data frame call ' thefiles' which will hold all the data at the end 
#right now just add the samecol into it as the first column 
thefiles <- data.frame(col=samecol, stringsAsFactors = FALSE)

#here is the loop to read all the files
#sorted it by the first column and remove the first column after sorted 
#combine the column to thefiles 
for (i in files) {
  the_file <- as.data.frame(read.csv(i, header = F))
  the_file <- the_file[order(the_file[[1]]),]
  the_file <- the_file[,-1]
  thefiles <- cbind(thefiles, the_file)
}

#right here we got all the data into one data frame, but not in the format we want
# cuz here we got unwanted column name here, if we transport it right now, we cannot remove it later 
#so remove the column name right now 
colnames(thefiles) <- NULL
#finally we can transport it rn 
thefiles <- t(thefiles)
#set the first row as the column name 
thefiles<- row_to_names(thefiles, row_number = 1)

dir.create("/Users/miasmacbook/Desktop/Group12/Document_made", showWarnings = FALSE)
write.csv(thefiles, "/Users/miasmacbook/Desktop/Group12/Document_made/merged_files.csv", row.names = FALSE)


dim(thefiles)

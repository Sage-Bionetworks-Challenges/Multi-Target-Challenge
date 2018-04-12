#! /usr/bin/env Rscript

library(synapser)
synLogin()



allCompounds = as.data.frame(synTableQuery('select * from syn12063758',includeRowIdAndRowVersion = FALSE))
problemInfo = as.data.frame(synTableQuery('select * from syn12081609',includeRowIdAndRowVersion = FALSE))

targets = union(grep(colnames(allCompounds),pattern = "CHEMB"), grep (colnames(allCompounds), pattern="pass"))
for (i in targets){ 
  toReplace = which(allCompounds[,i] == 'true')
  x = replace(allCompounds[,i],toReplace,TRUE)
  allCompounds[,i] = x
  }
for (i in targets){ 
  toReplace = which(allCompounds[,i] == 'false')
  x = replace(allCompounds[,i],toReplace,FALSE)
  allCompounds[,i] = as.logical(x)
}

prob1ReqBind = which(problemInfo$bind == 'true' & problemInfo$required == 'true' & problemInfo$problem == 1)
prob1AllBind = which(problemInfo$bind == 'true' & problemInfo$problem == 1)
prob2ReqBind = which(problemInfo$bind == 'true' & problemInfo$required == 'true' & problemInfo$problem == 2)
prob2AllBind = which(problemInfo$bind == 'true' & problemInfo$problem == 2)


# Does it pass for BINDING required for problem 1 (optional and anti-targets not considered)
colToCheck = which(colnames(allCompounds) %in% c(problemInfo$chembl_id[prob1ReqBind]))
rowProb1 = which(allCompounds$problem == 1)
passProbReq = allCompounds[which(allCompounds$problem == 1),colToCheck] == TRUE

passBindReq = rep(NA,nrow(allCompounds))
passBindReq[rowProb1] = passProbReq


# Does it pass for BINDING required for problem 2 (optional and anti-targets not considered)
colToCheck = which(colnames(allCompounds) %in% c(problemInfo$chembl_id[prob2ReqBind]))
rowProb2 = which(allCompounds$problem == 2)
passProbReq2 = allCompounds[which(allCompounds$problem == 2),colToCheck] == TRUE

# Check that will not overwrite any problem 1 values
intersect(rowProb1,rowProb2)
passBindReq[rowProb2] = passProbReq2[,1] & passProbReq2[,2]



# Does it pass for BINDING required and optional for problem 1 (anti-targets not considered)
colToCheck = which(colnames(allCompounds) %in% c(problemInfo$chembl_id[prob1AllBind]))
rowProb1 = which(allCompounds$problem == 1)
passProbAll = allCompounds[which(allCompounds$problem == 1),colToCheck] == TRUE

passBindAll = rep(NA,nrow(allCompounds))
passBindAll[rowProb1] = passProbAll[,1] & passProbAll[,2] & passProbAll[,3] & passProbAll[,4]


# Does it pass for BINDING Alluired and optional for problem 2 (anti-targets not considered)
colToCheck = which(colnames(allCompounds) %in% c(problemInfo$chembl_id[prob2AllBind]))
rowProb2 = which(allCompounds$problem == 2)
passProbAll2 = allCompounds[which(allCompounds$problem == 2),colToCheck] == TRUE

# Check that will not overwrite any problem 1 values
intersect(rowProb1,rowProb2)
passBindAll[rowProb2] = passProbAll2[,1] & passProbAll2[,2] & passProbAll2[,3] & passProbAll2[,4]


allCompounds$passBindReq = passBindReq
allCompounds$passBindAll = passBindAll



## How many compounds pass Tanimoto?

counts = matrix(NA,nrow = 4, ncol = 2)
colnames(counts) = c("pass", "fail")
rownames(counts) = c("required targets 1", "all targets 1", "required targets 2", "all targets 2")
sumCols = grep(colnames(allCompounds),pattern = "passBindReq")
counts[1,1] = sum(allCompounds[rowProb1,sumCols])
counts[3,1] = sum(allCompounds[rowProb2,sumCols])

sumCols = grep(colnames(allCompounds),pattern = "passBindAll")
counts[2,1] = sum(allCompounds[rowProb1,sumCols])
counts[4,1] = sum(allCompounds[rowProb2,sumCols])

# Count fail
counts[1,2] = length(rowProb1) - counts[1,1]
counts[3,2] = length(rowProb2) - counts[3,1]

counts[2,2] = length(rowProb1) - counts[2,1]
counts[4,2] = length(rowProb2) - counts[4,1]

op = par()
par(mar = c(10,4,4,2))
barplot(t(counts),beside = FALSE,legend.text = colnames(counts), las = 2, ylab = "number of submitted compounds", main = "Tanimoto score results")


## What is the average pass rate per target?
targets = grep(colnames(allCompounds),pattern = "CHEMB")
row1Pass = colSums(allCompounds[rowProb1,targets])
row2Pass = colSums(allCompounds[rowProb2,targets])

boxplot(row1Pass, row2Pass, ylab = "Number of compounds passing Tanimoto", main = "Solution uniqueness per target", names = c("problem 1", "problem 2"))



## How many compounds pass Lipinski

lipPSA = as.data.frame(synTableQuery('select submission_id, problem, solution, pass_lipinski, pass_polar_surface_area from syn12063027', includeRowIdAndRowVersion = FALSE))

targets = grep (colnames(lipPSA), pattern="pass")
for (i in targets){ 
  toReplace = which(lipPSA[,i] == 'true')
  x = replace(lipPSA[,i],toReplace,TRUE)
  lipPSA[,i] = x
}
for (i in targets){ 
  toReplace = which(lipPSA[,i] == 'false')
  x = replace(lipPSA[,i],toReplace,FALSE)
  lipPSA[,i] = as.logical(x)
}

counts = matrix(NA,nrow = 2, ncol = 2)
colnames(counts) = c("pass", "fail")
rownames(counts) = c("problem 1", "problem 2")
rowProb1 = which(lipPSA$problem == 1)
rowProb2 = which(lipPSA$problem == 2)
sumCols = grep(colnames(lipPSA),pattern = "lipinski")
counts[1,1] = sum(lipPSA[rowProb1,sumCols])
counts[2,1] = sum(lipPSA[rowProb2,sumCols])

# Count fail
counts[1,2] = length(rowProb1) - counts[1,1]
counts[2,2] = length(rowProb2) - counts[2,1]

op = par()
par(mar = c(10,4,4,2))
barplot(t(counts),beside = FALSE,legend.text = colnames(counts), las = 2, ylab = "number of submitted compounds", main = "Lipinski rule results")




## How many compounds pass PSA

counts = matrix(NA,nrow = 1, ncol = 2)
colnames(counts) = c("pass", "fail")
rownames(counts) = c("problem 2")
rowProb2 = which(lipPSA$problem == 2)
sumCols = grep(colnames(lipPSA),pattern = "polar")
counts[1,1] = sum(lipPSA[rowProb1,sumCols])

# Count fail
counts[1,2] = length(rowProb1) - counts[1,1]

op = par()
par(mar = c(10,4,4,2))
barplot(t(counts),beside = FALSE,legend.text = colnames(counts), las = 2, ylab = "number of submitted compounds", main = "PSA rule results")


## Total fraction passing filters

targets = union(union(grep(colnames(allCompounds),pattern = "CHEMB"), grep(colnames(allCompounds), pattern="for")), grep(colnames(allCompounds), pattern = "smiles"))
combinedData = merge(lipPSA, allCompounds[,-targets])

barplot(table(rowSums(combinedData[rowProb1,c(4,6)])), col = rainbow(5), main = "Success rate in passing criteria", ylab = "Number of compounds", xlab = "Number of criteria met")
barplot(table(rowSums(combinedData[rowProb2,c(4,5,6)])), col = rainbow(5), main = "Success rate in passing criteria", ylab = "Number of compounds", xlab = "Number of criteria met")

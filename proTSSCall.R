cArgs         <- commandArgs()
fileArg       <- which(startsWith(cArgs, "--file="))
argsArg       <- which(startsWith(cArgs, "--args"))
runningScript <- normalizePath(sapply(strsplit(cArgs[fileArg], "--file="), "[", 2))

tssPath       <- cArgs[argsArg+1]
readThreshold <- as.numeric(cArgs[argsArg+2])
tssList       <- read.table(tssPath, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

if(startsWith(runningScript, "/n/data1/cores/ntc")){
	scriptsPath <- "/n/data1/cores/ntc/scripts"
} else if(startsWith(runningScript, "/n/data1/hms/bcmp/adelman")){
	scriptsPath <- "/n/data1/hms/bcmp/adelman/Scripts/ntc"
} else{
	print("Warning: invoked script not on 'ntc' or 'adelman' paths")
	scriptsPath <- "/n/data1/cores/ntc/scripts"
}

setwd("bedGraphs")
system("mkdir -p forward reverse")
system("ln -rs *_3pr_forward.bedGraph forward")
system("ln -rs *_3pr_reverse.bedGraph reverse")

outPrefix <- sapply(strsplit(getwd(), "/"), tail, 2)[1]

setwd("forward")
system(
	paste(
		paste(scriptsPath, "AdelmanLab/NIH_scripts/bedgraphs2stdBedGraph/bedgraphs2stdBedGraph", sep = "/"),
		paste(outPrefix, "F", sep = "_")
	)
)

setwd("../reverse")
system(
	paste(
		paste(scriptsPath, "AdelmanLab/NIH_scripts/bedgraphs2stdBedGraph/bedgraphs2stdBedGraph", sep = "/"),
		paste(outPrefix, "R", sep = "_")
	)
)

setwd("../../")
tssPrefix <- sapply(
	strsplit(sapply(strsplit(tssPath, "/"), "[", length(strsplit(tssPath, "/")[[1]])), "_formakeheatmap.txt"),
	"[",
	1
)

system("mkdir pro_tss")

system(
	paste(
		scriptsPath,
		"/AdelmanLab/NIH_scripts/make_heatmap/make_heatmap -t 6 -l s -s s --nohead -p bedGraphs/forward/",
		outPrefix,
		"_F.bedGraph -m bedGraphs/reverse/",
		outPrefix,
		"_R.bedGraph -- ",
		tssPath,
		" pro_tss/",
		outPrefix,
		"_",
		tssPrefix,
		"_25mer_+-2kb.txt -2000 25 160", sep = ""
	)
)

tssMat <- read.table(
	paste("pro_tss/", outPrefix, "_", tssPrefix, "_25mer_+-2kb.txt", sep = ""),
	header    = TRUE,
	sep       = "\t",
	row.names = 1
)

# get TSS.+150 counts
tssReads    <- rowSums(tssMat[, which(colnames(tssMat) == "X0.24"):which(colnames(tssMat) == "X125.149")])
tssInactive <- names(tssReads)[tssReads <= readThreshold]

# print inactive TSSs
write.table(
	tssList[tssInactive, ],
	paste("pro_tss/", outPrefix, "_", tssPrefix, "_inactive_formakeheatmap.txt", sep = ""),
	quote     = FALSE,
	sep       = "\t",
	col.names = FALSE
)
tssActive <- names(tssReads)[tssReads > readThreshold]

# identify 1 dominant TSS per active gene
domTss <- vector("list")

for (i in tssActive) {
	if(!(tssList[i, "V7"] %in% names(domTss))) {
		domTss[[tssList[i, "V7"]]] <- i
	} else {
		# replace with current ENST, if more TSS-proximal reads

		if(tssReads[i] > tssReads[domTss[[tssList[i, "V7"]]]]) {
			domTss[[tssList[i, "V7"]]] <- i
		} else if (tssReads[i] == tssReads[domTss[[tssList[i, "V7"]]]]) {
			# if read tie, tiebreaker to more upstream ENST

			if(tssList[i, "V6"] == "+"){
				if(tssList[i, "V4"] < tssList[domTss[[tssList[i, "V7"]]], "V4"]) {
					domTss[[tssList[i, "V7"]]] <- i
				}
			} else if(tssList[i, "V6"] == "-"){
				if(tssList[i, "V5"] > tssList[domTss[[tssList[i, "V7"]]], "V5"]) {
					domTss[[tssList[i, "V7"]]] <- i
				}
			} else {
				stop("Unexpected strand symbol")
			}
		}
	} else {
	}
}

tssNonDom <- tssActive[!(tssActive %in% unlist(domTss))]
# print active, non-dominant TSSs
write.table(
	tssList[tssNonDom, ],
	paste("pro_tss/", outPrefix, "_", tssPrefix, "_non-dominant_formakeheatmap.txt", sep = ""),
	quote = FALSE,
	sep = "\t",
	col.names = FALSE
)

uniqueTss <- vector("list")
dupTss    <- vector("list")

# deduplicate TSSs
for (i in unlist(domTss)) {
	if(tssList[i, "V6"] == "+"){
		tssStart <- paste(tssList[i, "V3"], tssList[i, "V4"], "+", sep = ".")
	} else if(tssList[i, "V6"] == "-") {
		tssStart <- paste(tssList[i, "V3"], tssList[i, "V5"], "-", sep = ".")
	} else {
		stop("Unexpected strand symbol")
	}
	
	if(!(tssStart %in% names(uniqueTss))) {
		uniqueTss[[tssStart]] <- i
		dupTss[[tssStart]]    <- i
	} else {
		dupTss[[tssStart]] <- c(dupTss[[tssStart]], i)
		# if TSS exists, take longer transcript
		
		if(tssList[i, "V10"] > tssList[uniqueTss[[tssStart]], "V10"]){
			uniqueTss[[tssStart]] <- i
		} else if(tssList[i, "V10"] == tssList[uniqueTss[[tssStart]], "V10"]){
			# if TSS exists, and transcript lengths the same, take lower ENSG number
			
			digPos     <- regexpr("[[:digit:]]+", c(tssList[i, "V7"], tssList[uniqueTss[[tssStart]], "V7"]))
			geneDigits <- substr(c(tssList[i, "V7"], tssList[uniqueTss[[tssStart]], "V7"]), as.integer(digPos), as.integer(digPos) + (attributes(digPos)$match.length - 1))
			
			if(which.min(as.numeric(geneDigits)) == 1){
				uniqueTss[[tssStart]] <- i
			}
		} else {
		}
	}
}

# print dominant TSSs
write.table(
	tssList[rownames(tssList) %in% unlist(uniqueTss), ],
	paste0("pro_tss/", outPrefix, "_", tssPrefix, "_dominant_formakeheatmap.txt"),
	quote = FALSE,
	sep = "\t",
	col.names = FALSE
)

# print dominant gene bodies
system(
	paste(
		"awk -F'\t' -v OFS='\t' '{if($10>=400) {if($6==\"+\") print $1, $2, $3, $4+250, $5, $6, $7, $8, $9, $10; else print $1, $2, $3, $4, $5-250, $6, $7, $8, $9, $10}}'", 
		paste0("pro_tss/", outPrefix, "_", tssPrefix, "_dominant_formakeheatmap.txt"), ">", 
		paste0("pro_tss/", outPrefix, "_", tssPrefix, "_+250.TES.min400_formakeheatmap.txt")
	)
)

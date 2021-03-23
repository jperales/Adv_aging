
import os
import os.path
import re
import glob
import yaml

#--- Index
### Samples
SAMPLES = glob.glob("index/Samples/*.yaml")
SIDS = [re.sub("(AA|CA)_","", os.path.splitext(os.path.basename(k))[0]) for k in SAMPLES]

### Groups
GROUPS = glob.glob("index/Groups/*.yaml")
GIDS = [re.sub("RNA_","",os.path.splitext(os.path.basename(k))[0]) for k in GROUPS]

#--- Functions
def samplesFROMgroup(wildcards):
	yaml_fl = open("index/Groups/"+wildcards.region+"_"+wildcards.gid+".yaml")
	parsed_yaml = yaml.load(yaml_fl, Loader=yaml.FullLoader)
	samples = parsed_yaml["samples"].values()

	return samples

def RNA_getQCParams(wildcards):
	yaml_fl = open("index/Samples/"+wildcards.region+"_"+wildcards.id+".yaml")
	parsed_yaml = yaml.load(yaml_fl, Loader=yaml.FullLoader)
	QC_cutoff = parsed_yaml["QC_cutoff"]

	#NOTE: this will be deprecated when yaml template checks
	if ['MAX_percMT', 'MIN_nCells', 'MIN_nFeatures', 'MAX_nFeatures'] == list(QC_cutoff.keys()):
		QC_cutoff = list(QC_cutoff.values())
	else:
	 	sys.exit("ERROR: cutoff format is not correct")

	return QC_cutoff

def CPDBpval_FROMgroup(wildcards):
	yaml_fl = open("index/Groups/"+wildcards.gid+".yaml")
	parsed_yaml = yaml.load(yaml_fl, Loader=yaml.FullLoader)
	samples = parsed_yaml["samples"]
	fls = ["out/comm/Samples/"+k+"/"+wildcards.region+"/pvalues.txt" for k in samples]

	return fls

def CPDBmean_FROMgroup(wildcards):
	yaml_fl = open("index/Groups/"+wildcards.gid+".yaml")
	parsed_yaml = yaml.load(yaml_fl, Loader=yaml.FullLoader)
	samples = parsed_yaml["samples"]
	fls = ["out/comm/Samples/"+k+"/"+wildcards.region+"/means.txt" for k in samples]

	return fls

def CPDBsignif_FROMgroup(wildcards):
	yaml_fl = open("index/Groups/"+wildcards.gid+".yaml")
	parsed_yaml = yaml.load(yaml_fl, Loader=yaml.FullLoader)
	samples = parsed_yaml["samples"]
	fls = ["out/comm/Samples/"+k+"/"+wildcards.region+"/significant_means.txt" for k in samples]

	return fls





#--- Modules 
### All
rule all:
	input:
		expand("out/data2/Samples/{id}/{region}/{fl}", id=SIDS, region=["AA", "CA"], fl=["matrix.mtx", "barcodes.tsv", "features.tsv"]),
		expand("out/norm/Samples/{id}/{region}/logNormCounts{ext}", id=SIDS, region=["AA", "CA"], ext=["_matrix.mtx", "_barcodes.tsv", "_features.tsv"]),
		expand("out/norm/Samples/{id}/{region}/logNormSCT{ext}", id=SIDS, region=["AA", "CA"], ext=["_matrix.mtx", "_barcodes.tsv", "_features.tsv"]),
		expand("out/comm/Samples/{id}/{region}/cellphonedb_{suffix}", id=SIDS, region=["AA", "CA"], suffix=["meta.txt", "count.txt"]),
		expand("out/comm/Samples/{id}/{region}/{fl}", id=SIDS, region=["AA", "CA"], fl=["deconvoluted.txt", "means.txt", "pvalues.txt", "significant_means.txt"]),
		expand("out/comm/Groups/{gid}/{region}/merged.tsv", gid=GIDS, region=["AA", "CA"])

### Data process - hence data2
rule RNA_data_process:
	input:
		src = ["workflow/scripts/data_process.R"],
		fls = expand("data/GSE117715/{fl}", fl=["CellType_info.csv", "CellType_abbn.csv", "GSE117715_Cynomolgus_monkey_aging_artery_count.txt.gz"])
	output:
		mtx = expand("out/data2/Samples/{{id}}/{{region}}/{fl}", fl=["barcodes.tsv", "features.tsv", "matrix.mtx"]),
		rptQC = expand("out/QC/Samples/{{id}}/{{region}}/{fl}", fl=["barcodes_QC1.tsv", "features_QC1.tsv"])
	params:
		ID = lambda wildcards: wildcards.id,
		REG = lambda wildcards: wildcards.region,
		# Get QC params for cell filtering from sample YAML file
		QC_cutoff = RNA_getQCParams,
	message:
		"Data processing (filtering, report stats for QC)"
	conda:
		"workflow/envs/rna_process.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src[0]} {input.fls} {params.ID} {params.REG} {params.QC_cutoff} {output.mtx} {output.rptQC}"

### Normalization
rule RNA_norm_pooldeconv:
	input:
		src = ["workflow/scripts/norm_pooldeconv.R", "workflow/src/10Xmat.R"],
		mtx = expand("out/data2/Samples/{{id}}/{{region}}/{fl}", fl=["barcodes.tsv", "features.tsv", "matrix.mtx"]),
	output:
		mtx = expand("out/norm/Samples/{{id}}/{{region}}/logNormCounts_{suffix}", suffix=["barcodes.tsv", "features.tsv", "matrix.mtx"])
	message:
		"Normalizing by pool deconvolution (Lu et al 2016)"
	conda:
		"workflow/envs/osca.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src[0]} {input.mtx} {output.mtx}"

rule RNA_norm_sctransform:
	input:
		src = ["workflow/scripts/norm_sctransform.R", "workflow/src/10Xmat.R"],
		mtx = expand("out/data2/Samples/{{id}}/{{region}}/{fl}", fl=["barcodes.tsv", "features.tsv", "matrix.mtx"]),
	output:
		fls = expand("out/norm/Samples/{{id}}/{{region}}/logNormSCT_{suffix}", suffix=["barcodes.tsv", "features.tsv", "matrix.mtx", "scaled.mtx", "HGV.txt"])
	params:
		ID = lambda wildcards: wildcards.id
	message:
		"Normalizing using SCTransform (Hafemeister et al 2019)"
	conda:
		"workflow/envs/seurat.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src[0]} {input.mtx} {output.fls} {params.ID}"

### Communication
rule RNA_comm_prepare:
	input:
		src = ["workflow/scripts/comm_prepare.R", "workflow/src/10Xmat.R"],
		mtx = expand("out/norm/Samples/{{id}}/{{region}}/logNormSCT_{suffix}", suffix=["barcodes.tsv", "features.tsv", "matrix.mtx"])
	output:
		fls = expand("out/comm/Samples/{{id}}/{{region}}/cellphonedb_{suffix}", suffix=["meta.txt", "count.txt"])
	params:
		ID = lambda wildcards: wildcards.id,
		ANN = "Annotation.Level.1"
	conda:
		"workflow/envs/rna_process.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src[0]} {input.mtx} {params.ID} {params.ANN} {output.fls}"

#NOTE: implement cpdb params via yaml
rule RNA_comm_cpdb:
	input:
		fls = expand("out/comm/Samples/{{id}}/{{region}}/cellphonedb_{suffix}", suffix=["meta.txt", "count.txt"])
	output:
		expand("out/comm/Samples/{{id}}/{{region}}/{fl}", fl=["deconvoluted.txt", "means.txt", "pvalues.txt", "significant_means.txt"])
	params:
		outdir = lambda wildcards: "out/comm/Samples/"+wildcards.id+"/"+wildcards.region+"/"
	conda:
		"workflow/envs/cpdb.yaml"
	shell:
		"cellphonedb method statistical_analysis {input.fls} --iterations=100 --threshold=0.3 --counts-data=hgnc_symbol --output-path={params.outdir}"

rule RNA_comm_merge:
	input:
		src = "workflow/scripts/comm_merge.R",
		pval_fls = CPDBpval_FROMgroup,
		mean_fls = CPDBmean_FROMgroup,
		signif_fls = CPDBsignif_FROMgroup
	output:
		fl = "out/comm/Groups/{gid}/{region}/merged.tsv"
	params:
		pvals = lambda wildcards, input : ",".join(input.pval_fls),
		means = lambda wildcards, input : ",".join(input.mean_fls),
		signif = lambda wildcards, input : ",".join(input.signif_fls)
	conda:
		"workflow/envs/melt.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src} {params.pvals} {params.means} {params.signif} {output.fl}"
		

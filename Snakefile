
import os
import os.path
import re
import glob
import yaml
import pandas as pd
import numpy as np

#--- Index
### Samples
SAMPLES = glob.glob("index/Samples/*.yaml")
SIDS = [re.sub("(AA|CA)_","", os.path.splitext(os.path.basename(k))[0]) for k in SAMPLES]

### Groups
GROUPS = glob.glob("index/Groups/*.yaml")
GIDS = [re.sub("RNA_","",os.path.splitext(os.path.basename(k))[0]) for k in GROUPS]

### Contrasts
CONTRASTS = glob.glob("index/Contrasts/*.yaml")
CIDS = [os.path.splitext(os.path.basename(k))[0] for k in CONTRASTS]

#--- Functions

#- Index Group
def samplesFROMgroup(wildcards):
	yaml_fl = open("index/Groups/"+wildcards.gid+".yaml")
	parsed_yaml = yaml.load(yaml_fl, Loader=yaml.FullLoader)
	samples = parsed_yaml["samples"]

	return samples

#- Index Contrast
def groupsFROMcontrast(wildcards):
	groups = wildcards.cid.split("-")

	return groups

#- RNA processing
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

def sampleBARCODES_FROMgroup(wildcards):
	samples = samplesFROMgroup(wildcards)
	fls = ["out/data2/Samples/"+k+"/"+wildcards.region+"/barcodes.tsv" for k in samples]

	return fls

def sampleFEATURES_FROMgroup(wildcards):
	samples = samplesFROMgroup(wildcards)
	fls = ["out/data2/Samples/"+k+"/"+wildcards.region+"/features.tsv" for k in samples]

	return fls

def sampleMTX_FROMgroup(wildcards):
	samples = samplesFROMgroup(wildcards)
	fls = ["out/data2/Samples/"+k+"/"+wildcards.region+"/matrix.mtx" for k in samples]

	return fls

def groupBARCODES_FROMcontrast(wildcards):
	groups = groupsFROMcontrast(wildcards)
	fls = ["out/data2/Groups/"+k+"/"+wildcards.region+"/barcodes.tsv" for k in groups]

	return fls

def groupFEATURES_FROMcontrast(wildcards):
	groups = groupsFROMcontrast(wildcards)
	fls = ["out/data2/Groups/"+k+"/"+wildcards.region+"/features.tsv" for k in groups]

	return fls

def groupMTX_FROMcontrast(wildcards):
	groups = groupsFROMcontrast(wildcards)
	fls = ["out/data2/Groups/"+k+"/"+wildcards.region+"/matrix.mtx" for k in groups]

	return fls

def getAnn_FROMcontrast(wildcards):
	yaml_fl = open("index/Contrasts/"+wildcards.cid+".yaml")
	parsed_yaml = yaml.load(yaml_fl, Loader=yaml.FullLoader)
	ann = parsed_yaml["annotation"]

	return ann

def celltypes_FROMcontrast(wildcards):
	colData = pd.read_csv("out/data2/Contrasts/"+wildcards.cid+"/"+wildcards.region+"/barcodes.tsv", sep="\t")
	ann = getAnn_FROMcontrast(wildcards)
	clust = colData[ann].unique().tolist()
	return clust

#- CellPhoneDB
def CPDBpval_FROMgroup(wildcards):
	samples = samplesFROMgroup(wildcards)
	fls = ["out/comm/Samples/"+k+"/"+wildcards.region+"/pvalues.txt" for k in samples]

	return fls

def CPDBmean_FROMgroup(wildcards):
	samples = samplesFROMgroup(wildcards)
	fls = ["out/comm/Samples/"+k+"/"+wildcards.region+"/means.txt" for k in samples]

	return fls

def CPDBsignif_FROMgroup(wildcards):
	samples = samplesFROMgroup(wildcards)
	fls = ["out/comm/Samples/"+k+"/"+wildcards.region+"/significant_means.txt" for k in samples]

	return fls

def getINPUT_CrossTalker(wildcards):
	groups = groupsFROMcontrast(wildcards)
	fls = ["out/comm/Groups/"+k+"/"+wildcards.region+"/filtered_corrected.csv" for k in groups]
	return fls

# def pbMTX_FROMcontrast(CIDX):
# 	region = ["AA", "CA"]
# 	suffix = ["counts.tsv", "targets.tsv"]
# 	fls = []
# 	for k in CIDX:
# 		for j in region:
# 			dic1 = {"cid" : k, "region" : j}
# 			celltypes = celltypes_FROMcontrast(dic1)
# 			flsx = ["out/minibulk/"+k+"/"+j+"/"+cl for cl in celltypes]
# 			for fl in flsx:
# 				fls.append(fl)
# 	
# 	return fls




#--- Modules 
### All
rule all:
	input:
		expand("out/data2/Samples/{id}/{region}/{fl}", id=SIDS, region=["AA", "CA"], fl=["matrix.mtx", "barcodes.tsv", "features.tsv"]),
		expand("out/norm/Samples/{id}/{region}/logNormCounts{ext}", id=SIDS, region=["AA", "CA"], ext=["_matrix.mtx", "_barcodes.tsv", "_features.tsv"]),
		expand("out/norm/Samples/{id}/{region}/logNormSCT{ext}", id=SIDS, region=["AA", "CA"], ext=["_matrix.mtx", "_barcodes.tsv", "_features.tsv"]),
		expand("out/comm/Samples/{id}/{region}/cellphonedb_{suffix}", id=SIDS, region=["AA", "CA"], suffix=["meta.txt", "count.txt"]),
		expand("out/comm/Samples/{id}/{region}/{fl}", id=SIDS, region=["AA", "CA"], fl=["deconvoluted.txt", "means.txt", "pvalues.txt", "significant_means.txt"]),
		expand("out/comm/Groups/{gid}/{region}/{fl}", gid=GIDS, region=["AA", "CA"], fl=["merged.tsv", "legend.tsv"]),
		expand("out/comm/Groups/{gid}/{region}/combined.tsv", gid=GIDS, region=["AA", "CA"]),
		expand("out/comm/Groups/{gid}/{region}/{fl}", gid=GIDS, region=["AA", "CA"], fl=["combined_significant_means.tsv", "idx2cluster.tsv"]),
		expand("out/comm/Groups/{gid}/{region}/filtered_corrected.csv", gid=GIDS, region=["AA", "CA"]),
		expand("out/comm/Contrasts/{cid}/{region}/crosstalker/{prefix}_{cid}_{region}.html", cid=CIDS, region=["AA", "CA"], prefix=["Single", "Comparative"]),
		expand("out/comm/Contrasts/{cid}/{region}/crosstalker/{prefix}_{cid}_{region}.rds", cid=CIDS, region=["AA", "CA"], prefix=["data"]),
		expand("out/data2/Groups/{gid}/{region}/{fl}", gid=GIDS, region=["AA", "CA"], fl=["matrix.mtx", "barcodes.tsv", "features.tsv"]),
		expand("out/data2/Contrasts/{cid}/{region}/{fl}", cid=CIDS, region=["AA", "CA"], fl=["matrix.mtx", "barcodes.tsv", "features.tsv"]),
		dynamic(
			expand("out/minibulk/Contrasts/{cid}/{region}/{{cluster}}_{suffix}",
				cid=CIDS,
				region=["AA", "CA"],
				suffix=["counts.tsv", "samples.tsv"]
			)
		),
		dynamic(
			expand("out/minidge/Contrasts/{cid}/{region}/{{cluster}}_topTags.tsv",
				cid=CIDS,
				region=["AA", "CA"]
			)
		)

### Data process - hence data2
rule RNA_data2_process:
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

rule RNA_data2_merge_group:
	input:
		src = ["workflow/scripts/data_merge.R", "workflow/src/10Xmat.R"],
		barcodes_fls = sampleBARCODES_FROMgroup,
		features_fls = sampleFEATURES_FROMgroup,
		mtx_fls = sampleMTX_FROMgroup,
	output:
		mtx = expand("out/data2/Groups/{{gid}}/{{region}}/{fl}", fl=["barcodes.tsv", "features.tsv", "matrix.mtx"])
	params:
		barcodes = lambda wildcards, input : ",".join(input.barcodes_fls),
		features = lambda wildcards, input : ",".join(input.features_fls),
		mtx = lambda wildcards, input : ",".join(input.mtx_fls),
		GID = lambda wildcards: wildcards.gid,
		REG = lambda wildcards: wildcards.region
	conda:
		"workflow/envs/seurat.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src[0]} {params.barcodes} {params.features} {params.mtx} {params.GID} {params.REG} {output.mtx}"

rule RNA_data2_merge_contrast:
	input:
		src = ["workflow/scripts/data_merge.R", "workflow/src/10Xmat.R"],
		barcodes_fls = groupBARCODES_FROMcontrast,
		features_fls = groupFEATURES_FROMcontrast,
		mtx_fls = groupMTX_FROMcontrast,
	output:
		mtx = expand("out/data2/Contrasts/{{cid}}/{{region}}/{fl}", fl=["barcodes.tsv", "features.tsv", "matrix.mtx"])
	params:
		barcodes = lambda wildcards, input : ",".join(input.barcodes_fls),
		features = lambda wildcards, input : ",".join(input.features_fls),
		mtx = lambda wildcards, input : ",".join(input.mtx_fls),
		CID = lambda wildcards: wildcards.cid,
		REG = lambda wildcards: wildcards.region
	conda:
		"workflow/envs/seurat.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src[0]} {params.barcodes} {params.features} {params.mtx} {params.CID} {params.REG} {output.mtx}"

rule RNA_data2_pseudobulk_contrast:
	input:
		src = ["workflow/scripts/data_pseudobulk.R", "workflow/src/10Xmat.R"],
		mtx = expand("out/data2/Contrasts/{{cid}}/{{region}}/{fl}", fl=["barcodes.tsv", "features.tsv", "matrix.mtx"])
	output:
		dynamic("out/minibulk/Contrasts/{cid}/{region}/{cluster}_counts.tsv"),
		dynamic("out/minibulk/Contrasts/{cid}/{region}/{cluster}_samples.tsv")
	params:
		ann = getAnn_FROMcontrast,
		outdir = lambda wildcards: "out/minibulk/Contrasts/"+wildcards.cid+"/"+wildcards.region+"/"
	conda:
		"workflow/envs/osca.yaml"
	message:
		"Pseudobulking"
	shell:
		"Rscript --vanilla "
		"{input.src[0]} {input.mtx} {params.ann} {params.outdir}"

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
		fls = ["out/comm/Groups/{gid}/{region}/merged.tsv",
			"out/comm/Groups/{gid}/{region}/legend.tsv"]
	params:
		pvals = lambda wildcards, input : ",".join(input.pval_fls),
		means = lambda wildcards, input : ",".join(input.mean_fls),
		signif = lambda wildcards, input : ",".join(input.signif_fls)
	conda:
		"workflow/envs/melt.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src} {params.pvals} {params.means} {params.signif} {output.fls}"
		
rule RNA_comm_combine:
	input:
		src = "workflow/scripts/comm_combine.R",
		fl = "out/comm/Groups/{gid}/{region}/merged.tsv"
	output:
		fl = "out/comm/Groups/{gid}/{region}/combined.tsv"
#	params:
#		iterations = lambda wildcards, input : ",".join(input.pval_fls),
	conda:
		"workflow/envs/metapval.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src} {input.fl} {output.fl}"

rule RNA_comm_combineSignif:
	input:
		src = "workflow/scripts/comm_combineSignif.R",
		fls = ["out/comm/Groups/{gid}/{region}/combined.tsv",
			"out/comm/Groups/{gid}/{region}/legend.tsv"]
	output:
		fls = expand("out/comm/Groups/{{gid}}/{{region}}/{fl}", fl=["combined_significant_means.tsv", "idx2cluster.tsv"])
	conda:
		"workflow/envs/melt.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src} {input.fls} {output.fls}"

rule RNA_comm_cpdb2crosstalker:
	input:
		src = "workflow/scripts/comm_CPDB2CrossTalker.py",
		fls = ["out/comm/Groups/{gid}/{region}/combined_significant_means.tsv",
			"out/comm/Groups/{gid}/{region}/idx2cluster.tsv"]
	output:
		fl = "out/comm/Groups/{gid}/{region}/filtered_corrected.csv"
	conda:
		"workflow/envs/cpdb.yaml"
	shell:
		"python "
		"{input.src} {input.fls} {output.fl}"

rule RNA_comm_crosstalker:
	input:
		src = "workflow/scripts/comm_CrossTalker.R",
		fls = getINPUT_CrossTalker
	output:
		outdir = directory("out/comm/Contrasts/{cid}/{region}/crosstalker/"),
		reports = expand("out/comm/Contrasts/{{cid}}/{{region}}/crosstalker/{prefix}_{{cid}}_{{region}}.html", prefix=["Single", "Comparative"]),
		dat = "out/comm/Contrasts/{cid}/{region}/crosstalker/data_{cid}_{region}.rds"
	params:
		#suffix = lambda wildcards, output : re.sub("Single_","",os.path.splitext(os.path.basename(output.reports[0]))),
		suffix = lambda wildcards: wildcards.cid+"_"+wildcards.region,
		comparison = lambda wildcards: wildcards.cid,
		fls = lambda wildcards, input : ",".join(input.fls),
	shell:
		'set +eu '
		' && (test -d {output.outdir} || mkdir -p {output.outdir}) '
		' && . $(conda info --base)/etc/profile.d/conda.sh '
		' && conda activate envs/crosstalker '
		" && $CONDA_PREFIX/bin/Rscript --vanilla "
		"{input.src} {params.fls} {params.comparison} {params.suffix} {output.outdir}"

	
### pseudobulk testing
rule RNA_pseudobulk_dge:
	input:
		src = "workflow/scripts/dge_pseudobulk.R",
		cnt = "out/minibulk/Contrasts/{cid}/{region}/{cluster}_counts.tsv",
		trg = "out/minibulk/Contrasts/{cid}/{region}/{cluster}_samples.tsv"
	output:
		"out/minidge/Contrasts/{cid}/{region}/{cluster}_topTags.tsv",
	params:
		comparison = lambda wildcards: wildcards.cid
	conda:
		"workflow/envs/edger.yaml"
	message:
		"Pseudobulking"
	shell:
		"Rscript --vanilla "
		"{input.src} {input.cnt} {input.trg} {params.comparison} {output}"


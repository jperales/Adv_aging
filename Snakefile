
import os
import os.path
import re
import glob
import yaml
import pandas as pd
import numpy as np

#--- Gene Queries
#QUERIES = open("query.txt", 'r').read().splitlines()
QUERIES = ["PDGFRA"]


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

def RNA_getSeuratClustParams_FROMgroup(wildcards):
	yaml_fl = open("index/Groups/"+wildcards.gid+".yaml")
	parsed_yaml = yaml.load(yaml_fl, Loader=yaml.FullLoader)
	clust = ",".join(str(x) for x in parsed_yaml["res"])
	return clust 

def RNA_getSeuratClustParams_FROMcontrast(wildcards):
	yaml_fl = open("index/Contrasts/"+wildcards.cid+".yaml")
	parsed_yaml = yaml.load(yaml_fl, Loader=yaml.FullLoader)
	clust = ",".join(str(x) for x in parsed_yaml["res"])
	return clust 

def sampleUMIBARCODES_FROMgroup(wildcards):
	samples = samplesFROMgroup(wildcards)
	fls = ["out/data2/Samples/"+k+"/"+wildcards.region+"/barcodes.tsv" for k in samples]

	return fls

def sampleUMIFEATURES_FROMgroup(wildcards):
	samples = samplesFROMgroup(wildcards)
	fls = ["out/data2/Samples/"+k+"/"+wildcards.region+"/features.tsv" for k in samples]

	return fls

def sampleUMIMTX_FROMgroup(wildcards):
	samples = samplesFROMgroup(wildcards)
	fls = ["out/data2/Samples/"+k+"/"+wildcards.region+"/matrix.mtx" for k in samples]

	return fls

def groupUMIBARCODES_FROMcontrast(wildcards):
	groups = groupsFROMcontrast(wildcards)
	fls = ["out/data2/Groups/"+k+"/"+wildcards.region+"/barcodes.tsv" for k in groups]

	return fls

def groupUMIFEATURES_FROMcontrast(wildcards):
	groups = groupsFROMcontrast(wildcards)
	fls = ["out/data2/Groups/"+k+"/"+wildcards.region+"/features.tsv" for k in groups]

	return fls

def groupUMIMTX_FROMcontrast(wildcards):
	groups = groupsFROMcontrast(wildcards)
	fls = ["out/data2/Groups/"+k+"/"+wildcards.region+"/matrix.mtx" for k in groups]

	return fls

#-- Normalized data
def sampleNORMBARCODES_FROMgroup(wildcards):
	samples = samplesFROMgroup(wildcards)
	fls = ["out/norm/Samples/"+k+"/"+wildcards.region+"/"+wildcards.prefix+"_barcodes.tsv" for k in samples]

	return fls

def sampleNORMFEATURES_FROMgroup(wildcards):
	samples = samplesFROMgroup(wildcards)
	fls = ["out/norm/Samples/"+k+"/"+wildcards.region+"/"+wildcards.prefix+"_features.tsv" for k in samples]

	return fls

def sampleNORMMTX_FROMgroup(wildcards):
	samples = samplesFROMgroup(wildcards)
	fls = ["out/norm/Samples/"+k+"/"+wildcards.region+"/"+wildcards.prefix+"_matrix.mtx" for k in samples]

	return fls

def groupNORMBARCODES_FROMcontrast(wildcards):
	groups = groupsFROMcontrast(wildcards)
	fls = ["out/norm/Groups/"+k+"/"+wildcards.region+"/"+wildcards.prefix+"_barcodes.tsv" for k in groups]

	return fls

def groupNORMFEATURES_FROMcontrast(wildcards):
	groups = groupsFROMcontrast(wildcards)
	fls = ["out/norm/Groups/"+k+"/"+wildcards.region+"/"+wildcards.prefix+"_features.tsv" for k in groups]

	return fls

def groupNORMMTX_FROMcontrast(wildcards):
	groups = groupsFROMcontrast(wildcards)
	fls = ["out/norm/Groups/"+k+"/"+wildcards.region+"/"+wildcards.prefix+"_matrix.mtx" for k in groups]

	return fls

#-- Getting annotation
def getAnn_FROMcontrast(wildcards):
	yaml_fl = open("index/Contrasts/"+wildcards.cid+".yaml")
	parsed_yaml = yaml.load(yaml_fl, Loader=yaml.FullLoader)
	ann = parsed_yaml["annotation"]

	return ann

def getAnnFile_FROMsample(wildcards):
	yaml_fl = open("index/Samples/"+wildcards.region+"_"+wildcards.id+".yaml")
	parsed_yaml = yaml.load(yaml_fl, Loader=yaml.FullLoader)
	gr = parsed_yaml["biotype"]

	fl = "out/ann/Groups/"+gr+"/"+wildcards.region+"/"+"logNormSCT_harmony_barcodes.tsv"

	return fl


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

def getINPUT2_CrossTalker(wildcards):
	groups = groupsFROMcontrast(wildcards)
	fls = ["out/cpdb/Groups/"+k+"/"+wildcards.region+"/filtered_corrected.csv" for k in groups]
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
		expand("out/data2/Groups/{gid}/{region}/{fl}", gid=GIDS, region=["AA", "CA"], fl=["matrix.mtx", "barcodes.tsv", "features.tsv"]),
		expand("out/data2/Contrasts/{cid}/{region}/{fl}", cid=CIDS, region=["AA", "CA"], fl=["matrix.mtx", "barcodes.tsv", "features.tsv"]),
		expand("out/norm/Samples/{id}/{region}/logNormCounts{ext}", id=SIDS, region=["AA", "CA"], ext=["_matrix.mtx", "_barcodes.tsv", "_features.tsv"]),
		expand("out/norm/Samples/{id}/{region}/logNormSCT{ext}", id=SIDS, region=["AA", "CA"], ext=["_matrix.mtx", "_barcodes.tsv", "_features.tsv"]),
		expand("out/norm/Groups/{gid}/{region}/logNormSCT_{fl}", 
				gid=GIDS, region=["AA", "CA"],
				fl=["barcodes.tsv", "features.tsv", "counts.mtx", "data.mtx","scaled.tsv", "HGV.txt"]),
		expand("out/norm/Contrasts/{cid}/{region}/logNormSCT_{fl}", 
				cid=CIDS, region=["AA", "CA"],
				fl=["barcodes.tsv", "features.tsv", "counts.mtx", "data.mtx","scaled.tsv", "HGV.txt"]),
		expand("out/dim/Groups/{gid}/{region}/logNormSCT_PCA_{suffix}", gid=GIDS, region=["AA", "CA"], suffix=["embeddings.tsv", "loadings.tsv", "projected.tsv", "stdev.txt", "nPCs.txt", "horn.pdf"]),
		expand("out/dim/Contrasts/{cid}/{region}/logNormSCT_PCA_{suffix}", cid=CIDS, region=["AA", "CA"], suffix=["embeddings.tsv", "loadings.tsv", "projected.tsv", "stdev.txt", "nPCs.txt", "horn.pdf"]),
		expand("out/dim/Groups/{gid}/{region}/logNormSCT_harmony_{suffix}", gid=GIDS, region=["AA", "CA"], suffix=["embeddings.tsv", "loadings.tsv", "projected.tsv", "stdev.txt", "nPCs.txt"]),
		expand("out/dim/Contrasts/{cid}/{region}/logNormSCT_harmony_{suffix}", cid=CIDS, region=["AA", "CA"], suffix=["embeddings.tsv", "loadings.tsv", "projected.tsv", "stdev.txt", "nPCs.txt"]),
		expand("out/clust/Groups/{gid}/{region}/logNormSCT_harmony_{suffix}", gid=GIDS, region=["AA", "CA"], suffix=["GraphNN.rds", "GraphSNN.rds", "Idents.tsv"]),
		expand("out/clust/Contrasts/{cid}/{region}/logNormSCT_harmony_{suffix}", cid=CIDS, region=["AA", "CA"], suffix=["GraphNN.rds", "GraphSNN.rds", "Idents.tsv"]),
#		expand("out/ann/Samples/{id}/{region}/logNormSCT_harmony_{suffix}", id=SIDS, region=["AA", "CA"], suffix=["barcodes.tsv"]),
		expand("out/ann/Groups/{gid}/{region}/logNormSCT_harmony_{suffix}", gid=GIDS, region=["AA", "CA"], suffix=["barcodes.tsv", "vis.pdf"]),
		expand("out/ann/Contrasts/{cid}/{region}/logNormSCT_harmony_{suffix}", cid=CIDS, region=["AA", "CA"], suffix=["barcodes.tsv", "vis.pdf"]),
#		expand("out/comm/Samples/{id}/{region}/cellphonedb_{suffix}", id=SIDS, region=["AA", "CA"], suffix=["meta.txt", "count.txt"]),
#		expand("out/comm/Samples/{id}/{region}/{fl}", id=SIDS, region=["AA", "CA"], fl=["deconvoluted.txt", "means.txt", "pvalues.txt", "significant_means.txt"]),
#		expand("out/comm/Groups/{gid}/{region}/{fl}", gid=GIDS, region=["AA", "CA"], fl=["merged.tsv", "legend.tsv"]),
#		expand("out/comm/Groups/{gid}/{region}/combined.tsv", gid=GIDS, region=["AA", "CA"]),
#		expand("out/comm/Groups/{gid}/{region}/{fl}", gid=GIDS, region=["AA", "CA"], fl=["combined_significant_means.tsv", "idx2cluster.tsv"]),
#		expand("out/comm/Groups/{gid}/{region}/filtered_corrected.csv", gid=GIDS, region=["AA", "CA"]),
#		expand("out/comm/Contrasts/{cid}/{region}/crosstalker/{prefix}_{cid}_{region}.html", cid=CIDS, region=["AA", "CA"], prefix=["Single", "Comparative"]),
#		expand("out/comm/Contrasts/{cid}/{region}/crosstalker/{prefix}_{cid}_{region}.rds", cid=CIDS, region=["AA", "CA"], prefix=["data"]),

###		expand("out/cpdb/Samples/{id}/{region}/cellphonedb_{suffix}", id=SIDS, region=["AA", "CA"], suffix=["meta.txt", "count.txt"]),
#		expand("out/cpdb/Groups/{gid}/{region}/cellphonedb_{suffix}", gid=GIDS, region=["AA", "CA"], suffix=["meta.txt", "count.txt"]),
###		expand("out/cpdb/Samples/{id}/{region}/{fl}", id=SIDS, region=["AA", "CA"], fl=["deconvoluted.txt", "means.txt", "pvalues.txt", "significant_means.txt"]),
#		expand("out/cpdb/Groups/{gid}/{region}/{fl}", gid=GIDS, region=["AA", "CA"], fl=["deconvoluted.txt", "means.txt", "pvalues.txt", "significant_means.txt"]),
#		expand("out/cpdb/Groups/{gid}/{region}/{fl}", gid=GIDS, region=["AA", "CA"], fl=["significant_means.tsv", "idx2cluster.tsv"]),
#		expand("out/cpdb/Groups/{gid}/{region}/filtered_corrected.csv", gid=GIDS, region=["AA", "CA"]),
#		expand("out/cpdb/Contrasts/{cid}/{region}/crosstalker/{prefix}_{cid}_{region}.html", cid=CIDS, region=["AA", "CA"], prefix=["Single", "Comparative"]),
#		expand("out/cpdb/Contrasts/{cid}/{region}/crosstalker/{prefix}_{cid}_{region}.rds", cid=CIDS, region=["AA", "CA"], prefix=["data"]),
		expand("out/report/{cid}/{query}/violin_{cid}_{region}.{ext}", cid=CIDS, region=["AA", "CA"], ext=["pdf", "png", "tiff"], query=QUERIES),
		expand("out/export/SeuratObject/Zhang2020_AgingArteryAtlasCynomolgusMonkey_{cid}_{region}.rds", cid=CIDS, region=["AA", "CA"])
#		dynamic(
#			expand("out/minibulk/Contrasts/{cid}/{region}/{{cluster}}_{suffix}",
#				cid=CIDS,
#				region=["AA", "CA"],
#				suffix=["counts.tsv", "samples.tsv"]
#			)
#		),
#		dynamic(
#			expand("out/minidge/Contrasts/{cid}/{region}/{{cluster}}_topTags.tsv",
#				cid=CIDS,
#				region=["AA", "CA"]
#			)
#		),

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
		barcodes_fls = sampleUMIBARCODES_FROMgroup,
		features_fls = sampleUMIFEATURES_FROMgroup,
		mtx_fls = sampleUMIMTX_FROMgroup,
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
		barcodes_fls = groupUMIBARCODES_FROMcontrast,
		features_fls = groupUMIFEATURES_FROMcontrast,
		mtx_fls = groupUMIMTX_FROMcontrast,
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

### Dimensionality reduction shorten as dim
rule RNA_dim_pca_group:
	input:
		src = ["workflow/scripts/dim_pca_merged.R", "workflow/src/10Xmat.R"],
		umi = expand("out/data2/Groups/{{gid}}/{{region}}/{fl}", fl=["barcodes.tsv", "features.tsv", "matrix.mtx"]),
		sct = expand("out/norm/Groups/{{gid}}/{{region}}/logNormSCT_{fl}", 
				fl=["barcodes.tsv", "features.tsv", "counts.mtx", "data.mtx", "scaled.tsv", "HGV.txt"])
	output:
		fls = expand("out/dim/Groups/{{gid}}/{{region}}/logNormSCT_PCA_{suffix}", suffix=["embeddings.tsv", "loadings.tsv", "projected.tsv", "stdev.txt", "nPCs.txt", "horn.pdf"])
	params:
		GID = lambda wildcards: wildcards.gid,
		REG = lambda wildcards: wildcards.region
	conda:
		"workflow/envs/seurat.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src[0]} {input.umi} {input.sct} {params.GID} {params.REG} {output.fls}"

rule RNA_dim_pca_contrast:
	input:
		src = ["workflow/scripts/dim_pca_merged.R", "workflow/src/10Xmat.R"],
		umi = expand("out/data2/Contrasts/{{cid}}/{{region}}/{fl}", fl=["barcodes.tsv", "features.tsv", "matrix.mtx"]),
		sct = expand("out/norm/Contrasts/{{cid}}/{{region}}/logNormSCT_{fl}", 
				fl=["barcodes.tsv", "features.tsv", "counts.mtx", "data.mtx", "scaled.tsv", "HGV.txt"])
	output:
		fls = expand("out/dim/Contrasts/{{cid}}/{{region}}/logNormSCT_PCA_{suffix}", suffix=["embeddings.tsv", "loadings.tsv", "projected.tsv", "stdev.txt", "nPCs.txt", "horn.pdf"])
	params:
		GID = lambda wildcards: wildcards.cid,
		REG = lambda wildcards: wildcards.region
	conda:
		"workflow/envs/seurat.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src[0]} {input.umi} {input.sct} {params.GID} {params.REG} {output.fls}"

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

rule RNA_norm_merge_group:
	input:
		src = ["workflow/scripts/norm_merge.R", "workflow/src/10Xmat.R"],
		mtx = expand("out/data2/Groups/{{gid}}/{{region}}/{fl}", fl=["barcodes.tsv", "features.tsv", "matrix.mtx"])
		#UMIbarcodes_fls = sampleUMIBARCODES_FROMgroup,
		#UMIfeatures_fls = sampleUMIFEATURES_FROMgroup,
		#UMImtx_fls = sampleUMIMTX_FROMgroup,
		#NORMbarcodes_fls = sampleNORMBARCODES_FROMgroup,
		#NORMfeatures_fls = sampleNORMFEATURES_FROMgroup,
		#NORMmtx_fls = sampleNORMMTX_FROMgroup,
	output:
		mtx = expand("out/norm/Groups/{{gid}}/{{region}}/logNormSCT_{fl}", fl=["barcodes.tsv", "features.tsv", "counts.mtx", "data.mtx", "scaled.tsv", "HGV.txt"])
	params:
		#		UMIbarcodes = lambda wildcards, input : ",".join(input.UMIbarcodes_fls),
		#		UMIfeatures = lambda wildcards, input : ",".join(input.UMIfeatures_fls),
		#		UMImtx = lambda wildcards, input : ",".join(input.UMImtx_fls),
		#		NORMbarcodes = lambda wildcards, input : ",".join(input.NORMbarcodes_fls),
		#		NORMfeatures = lambda wildcards, input : ",".join(input.NORMfeatures_fls),
		#		NORMmtx = lambda wildcards, input : ",".join(input.NORMmtx_fls),
		GID = lambda wildcards: wildcards.xid,
		REG = lambda wildcards: wildcards.region
	conda:
		"workflow/envs/seurat.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src[0]} "
		"{input.mtx} "
		#"{params.UMIbarcodes} {params.UMIfeatures} {params.UMImtx} "
		#"{params.NORMbarcodes} {params.NORMfeatures} {params.NORMmtx} "
		"{params.GID} {params.REG} {output.mtx}"

rule RNA_norm_merge_contrast:
	input:
		src = ["workflow/scripts/norm_merge.R", "workflow/src/10Xmat.R"],
		mtx = expand("out/data2/Contrasts/{{cid}}/{{region}}/{fl}", fl=["barcodes.tsv", "features.tsv", "matrix.mtx"])
	output:
		mtx = expand("out/norm/Contrasts/{{cid}}/{{region}}/logNormSCT_{fl}", fl=["barcodes.tsv", "features.tsv", "counts.mtx", "data.mtx", "scaled.tsv", "HGV.txt"])
	params:
		GID = lambda wildcards: wildcards.cid,
		REG = lambda wildcards: wildcards.region
	conda:
		"workflow/envs/seurat.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src[0]} "
		"{input.mtx} "
		"{params.GID} {params.REG} {output.mtx}"

### Integration
rule RNA_int_harmony_group:
	input:
		src = ["workflow/scripts/int_harmony_group.R", "workflow/src/10Xmat.R"],
		umi = expand("out/data2/Groups/{{gid}}/{{region}}/{fl}", fl=["barcodes.tsv", "features.tsv", "matrix.mtx"]),
		sct = expand("out/norm/Groups/{{gid}}/{{region}}/logNormSCT_{fl}", 
				fl=["barcodes.tsv", "features.tsv", "counts.mtx", "data.mtx", "scaled.tsv", "HGV.txt"]),
		pca = expand("out/dim/Groups/{{gid}}/{{region}}/logNormSCT_PCA_{suffix}", suffix=["embeddings.tsv", "loadings.tsv", "projected.tsv", "stdev.txt", "nPCs.txt"])
	output:
		fls = expand("out/dim/Groups/{{gid}}/{{region}}/logNormSCT_harmony_{suffix}", suffix=["embeddings.tsv", "loadings.tsv", "projected.tsv", "stdev.txt", "nPCs.txt"])
	params:
		GID = lambda wildcards: wildcards.gid,
		REG = lambda wildcards: wildcards.region
	shell:
		'set +eu '
		' && . $(conda info --base)/etc/profile.d/conda.sh '
		' && conda activate envs/harmony '
		" && $CONDA_PREFIX/bin/Rscript --vanilla "
		"{input.src[0]} {input.umi} {input.sct} {input.pca} {params.GID} {params.REG} {output.fls}"

rule RNA_int_harmony_contrast:
	input:
		#NOTE: Shall we correct for group in harmony? Hence individual script for it. See prev rule src
		src = ["workflow/scripts/int_harmony_group.R", "workflow/src/10Xmat.R"],
		umi = expand("out/data2/Contrasts/{{cid}}/{{region}}/{fl}", fl=["barcodes.tsv", "features.tsv", "matrix.mtx"]),
		sct = expand("out/norm/Contrasts/{{cid}}/{{region}}/logNormSCT_{fl}", 
				fl=["barcodes.tsv", "features.tsv", "counts.mtx", "data.mtx", "scaled.tsv", "HGV.txt"]),
		pca = expand("out/dim/Contrasts/{{cid}}/{{region}}/logNormSCT_PCA_{suffix}", suffix=["embeddings.tsv", "loadings.tsv", "projected.tsv", "stdev.txt", "nPCs.txt"])
	output:
		fls = expand("out/dim/Contrasts/{{cid}}/{{region}}/logNormSCT_harmony_{suffix}", suffix=["embeddings.tsv", "loadings.tsv", "projected.tsv", "stdev.txt", "nPCs.txt"])
	params:
		GID = lambda wildcards: wildcards.cid,
		REG = lambda wildcards: wildcards.region
	shell:
		'set +eu '
		' && . $(conda info --base)/etc/profile.d/conda.sh '
		' && conda activate envs/harmony '
		" && $CONDA_PREFIX/bin/Rscript --vanilla "
		"{input.src[0]} {input.umi} {input.sct} {input.pca} {params.GID} {params.REG} {output.fls}"

### Clustering
rule RNA_clust_harmony_group:
	input:
		src = ["workflow/scripts/clust_SNN.R", "workflow/src/10Xmat.R"],
		umi = expand("out/data2/Groups/{{gid}}/{{region}}/{fl}", fl=["barcodes.tsv", "features.tsv", "matrix.mtx"]),
		sct = expand("out/norm/Groups/{{gid}}/{{region}}/logNormSCT_{fl}", 
				fl=["barcodes.tsv", "features.tsv", "counts.mtx", "data.mtx", "scaled.tsv", "HGV.txt"]),
		harmony = expand("out/dim/Groups/{{gid}}/{{region}}/logNormSCT_harmony_{suffix}", suffix=["embeddings.tsv", "loadings.tsv", "projected.tsv", "stdev.txt", "nPCs.txt"])
	output:
		fls = expand("out/clust/Groups/{{gid}}/{{region}}/logNormSCT_harmony_{suffix}", suffix=["GraphNN.rds", "GraphSNN.rds", "Idents.tsv"])
	params:
		GID = lambda wildcards: wildcards.gid,
		REG = lambda wildcards: wildcards.region,
		res = RNA_getSeuratClustParams_FROMgroup,
		red = "harmony"
	conda:
		"workflow/envs/seurat.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src[0]} {input.umi} {input.sct} {input.harmony} {params.GID} {params.REG} {params.res} {params.red} {output.fls}"

#NOTE: depcreated: annotation is propagated from integrated 2-group (contrasts)
##rule RNA_clust_subclust_group:
##	input:
##		src = ["workflow/scripts/clust_subclust.R", "workflow/src/10Xmat.R"],
##		umi = expand("out/data2/Groups/{{gid}}/{{region}}/{fl}", fl=["barcodes.tsv", "features.tsv", "matrix.mtx"]),
##		sct = expand("out/norm/Groups/{{gid}}/{{region}}/logNormSCT_{fl}", 
##				fl=["barcodes.tsv", "features.tsv", "counts.mtx", "data.mtx", "scaled.tsv", "HGV.txt"]),
##		harmony = expand("out/dim/Groups/{{gid}}/{{region}}/logNormSCT_harmony_{suffix}", suffix=["embeddings.tsv", "loadings.tsv", "projected.tsv", "stdev.txt", "nPCs.txt"]),
##		graph = expand("out/clust/Groups/{{gid}}/{{region}}/logNormSCT_harmony_{suffix}", suffix=["GraphNN.rds", "GraphSNN.rds"])
##	output:
##		fls = expand("out/ann/Groups/{{gid}}/{{region}}/logNormSCT_harmony_{suffix}", suffix=["barcodes.tsv", "vis.pdf"]),
##	params:
##		GID = lambda wildcards: wildcards.gid,
##		REG = lambda wildcards: wildcards.region,
##		#res = RNA_getSeuratClustParams_FROMgroup,
##		res = "0.5",
##		red = "harmony"
##	conda:
##		"workflow/envs/seurat.yaml"
##	shell:
##		"Rscript --vanilla "
##		"{input.src[0]} {input.umi} {input.sct} {input.harmony} {input.graph} {params.GID} {params.REG} {params.res} {params.red} {output.fls}"

rule RNA_clust_harmony_contrast:
	input:
		src = ["workflow/scripts/clust_SNN.R", "workflow/src/10Xmat.R"],
		umi = expand("out/data2/Contrasts/{{cid}}/{{region}}/{fl}", fl=["barcodes.tsv", "features.tsv", "matrix.mtx"]),
		sct = expand("out/norm/Contrasts/{{cid}}/{{region}}/logNormSCT_{fl}", 
				fl=["barcodes.tsv", "features.tsv", "counts.mtx", "data.mtx", "scaled.tsv", "HGV.txt"]),
		harmony = expand("out/dim/Contrasts/{{cid}}/{{region}}/logNormSCT_harmony_{suffix}", suffix=["embeddings.tsv", "loadings.tsv", "projected.tsv", "stdev.txt", "nPCs.txt"])
	output:
		fls = expand("out/clust/Contrasts/{{cid}}/{{region}}/logNormSCT_harmony_{suffix}", suffix=["GraphNN.rds", "GraphSNN.rds", "Idents.tsv"])
	params:
		GID = lambda wildcards: wildcards.cid,
		REG = lambda wildcards: wildcards.region,
		res = RNA_getSeuratClustParams_FROMgroup,
		red = "harmony"
	conda:
		"workflow/envs/seurat.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src[0]} {input.umi} {input.sct} {input.harmony} {params.GID} {params.REG} {params.res} {params.red} {output.fls}"

rule RNA_clust_subclust_contrast:
	input:
		src = ["workflow/scripts/clust_subclust.R", "workflow/src/10Xmat.R"],
		umi = expand("out/data2/Contrasts/{{cid}}/{{region}}/{fl}", fl=["barcodes.tsv", "features.tsv", "matrix.mtx"]),
		sct = expand("out/norm/Contrasts/{{cid}}/{{region}}/logNormSCT_{fl}", 
				fl=["barcodes.tsv", "features.tsv", "counts.mtx", "data.mtx", "scaled.tsv", "HGV.txt"]),
		harmony = expand("out/dim/Contrasts/{{cid}}/{{region}}/logNormSCT_harmony_{suffix}", suffix=["embeddings.tsv", "loadings.tsv", "projected.tsv", "stdev.txt", "nPCs.txt"]),
		graph = expand("out/clust/Contrasts/{{cid}}/{{region}}/logNormSCT_harmony_{suffix}", suffix=["GraphNN.rds", "GraphSNN.rds"])
	output:
		fls = expand("out/ann/Contrasts/{{cid}}/{{region}}/logNormSCT_harmony_{suffix}", suffix=["barcodes.tsv", "vis.pdf"]),
	params:
		GID = lambda wildcards: wildcards.cid,
		REG = lambda wildcards: wildcards.region,
		#res = RNA_getSeuratClustParams_FROMcontrast,
		res = "0.1",
		red = "harmony"
	conda:
		"workflow/envs/seurat.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src[0]} {input.umi} {input.sct} {input.harmony} {input.graph} {params.GID} {params.REG} {params.res} {params.red} {output.fls}"


rule RNA_ann_sample:
	input:
		src = ["workflow/scripts/ann_sample.R", "workflow/src/10Xmat.R"],
		norm = "out/norm/Samples/{id}/{region}/logNormSCT_barcodes.tsv", 
#		ann = getAnnFile_FROMsample
		ann = "out/ann/Contrasts/old-young/{region}/logNormSCT_harmony_barcodes.tsv"
	output:
		mtx_1 = "out/ann/Samples/{id}/{region}/logNormSCT_harmony_barcodes.tsv",
	conda:
		"workflow/envs/rna_process.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src[0]} {input.norm} {input.ann} {output.mtx_1}"

rule RNA_ann_group:
	input:
		src = ["workflow/scripts/ann_sample.R", "workflow/src/10Xmat.R"],
		norm = "out/norm/Groups/{gid}/{region}/logNormSCT_barcodes.tsv", 
#		ann = getAnnFile_FROMsample
		ann = "out/ann/Contrasts/old-young/{region}/logNormSCT_harmony_barcodes.tsv"
	output:
		mtx_1 = "out/ann/Groups/{gid}/{region}/logNormSCT_harmony_barcodes.tsv",
	conda:
		"workflow/envs/rna_process.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src[0]} {input.norm} {input.ann} {output.mtx_1}"

### Communication
#NOTE: Communication module is composed of two branches: 
# - comm: classical cellphonedb runs and
# - cpdb: single-sample cellphonedb runs followed by meta-analysis. 
# Both finish with crosstalker # for a differential communication between 2-group.
rule RNA_comm_prepare:
	input:
		src = ["workflow/scripts/comm_prepare.R", "workflow/src/10Xmat.R"],
		#mtx = expand("out/norm/Samples/{{id}}/{{region}}/logNormSCT_{suffix}", suffix=["barcodes.tsv", "features.tsv", "matrix.mtx"])
		mtx_1 = "out/ann/Samples/{id}/{region}/logNormSCT_harmony_barcodes.tsv",
		mtx_23 = expand("out/norm/Samples/{{id}}/{{region}}/logNormSCT_{suffix}", suffix=["features.tsv", "matrix.mtx"])
	output:
		fls = expand("out/comm/Samples/{{id}}/{{region}}/cellphonedb_{suffix}", suffix=["meta.txt", "count.txt"])
	params:
		ID = lambda wildcards: wildcards.id,
		ANN = "Annotation.Level.2"
	conda:
		"workflow/envs/rna_process.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src[0]} {input.mtx_1} {input.mtx_23} {params.ID} {params.ANN} {output.fls}"

rule RNA_cpdb_prepare:
	input:
		src = ["workflow/scripts/comm_prepare.R", "workflow/src/10Xmat.R"],
		#mtx = expand("out/norm/Samples/{{id}}/{{region}}/logNormSCT_{suffix}", suffix=["barcodes.tsv", "features.tsv", "matrix.mtx"])
		mtx_1 = "out/ann/Groups/{gid}/{region}/logNormSCT_harmony_barcodes.tsv",
		mtx_23 = expand("out/norm/Groups/{{gid}}/{{region}}/logNormSCT_{suffix}", suffix=["features.tsv", "data.mtx"])
	output:
		fls = expand("out/cpdb/Groups/{{gid}}/{{region}}/cellphonedb_{suffix}", suffix=["meta.txt", "count.txt"])
	params:
		ID = lambda wildcards: wildcards.gid,
		ANN = "Annotation.Level.2"
	conda:
		"workflow/envs/rna_process.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src[0]} {input.mtx_1} {input.mtx_23} {params.ID} {params.ANN} {output.fls}"

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

rule RNA_cpdb_cpdb:
	input:
		fls = expand("out/cpdb/Groups/{{gid}}/{{region}}/cellphonedb_{suffix}", suffix=["meta.txt", "count.txt"])
	output:
		expand("out/cpdb/Groups/{{gid}}/{{region}}/{fl}", fl=["deconvoluted.txt", "means.txt", "pvalues.txt", "significant_means.txt"])
	params:
		outdir = lambda wildcards: "out/cpdb/Groups/"+wildcards.gid+"/"+wildcards.region+"/"
	conda:
		"workflow/envs/cpdb.yaml"
	shell:
		"cellphonedb method statistical_analysis {input.fls} --iterations=1000 --threshold=0.3 --counts-data=hgnc_symbol --output-path={params.outdir}"

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

#NOTE: eq to RNA_comm_combineSignif
rule RNA_cpdb_idx2cluster:
	input:
		src = "workflow/scripts/comm_getIdx2cluster.R",
		fl = "out/cpdb/Groups/{gid}/{region}/significant_means.txt"
	output:
		fls = ["out/cpdb/Groups/{gid}/{region}/significant_means.tsv",
			"out/cpdb/Groups/{gid}/{region}/idx2cluster.tsv"]
	conda:
		"workflow/envs/rna_process.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src} {input.fl} {output.fls}"

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

rule RNA_cpdb_cpdb2crosstalker:
	input:
		src = "workflow/scripts/comm_CPDB2CrossTalker.py",
		fls = ["out/cpdb/Groups/{gid}/{region}/significant_means.tsv",
			"out/cpdb/Groups/{gid}/{region}/idx2cluster.tsv"]
	output:
		fl = "out/cpdb/Groups/{gid}/{region}/filtered_corrected.csv"
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

rule RNA_cpdb_crosstalker:
	input:
		src = "workflow/scripts/comm_CrossTalker.R",
		fls = getINPUT2_CrossTalker
	output:
		outdir = directory("out/cpdb/Contrasts/{cid}/{region}/crosstalker/"),
		reports = expand("out/cpdb/Contrasts/{{cid}}/{{region}}/crosstalker/{prefix}_{{cid}}_{{region}}.html", prefix=["Single", "Comparative"]),
		dat = "out/cpdb/Contrasts/{cid}/{region}/crosstalker/data_{cid}_{region}.rds"
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

### Visualization
rule vis_vlnplot:
	input:
		src = ["workflow/scripts/vis_violin.R", "workflow/src/10Xmat.R"],
		barcodes = "out/ann/Contrasts/{cid}/{region}/logNormSCT_harmony_barcodes.tsv",
		norm_mtx = expand("out/norm/Contrasts/{{cid}}/{{region}}/logNormSCT_{suffix}", suffix=["features.tsv", "counts.mtx", "data.mtx"])
	output:
		"out/report/{cid}/{query}/violin_{cid}_{region}.{ext}"
	params:
		ds = lambda wildcards: wildcards.cid,
		query = lambda wildcards: wildcards.query,
		ann = "Annotation.Level.2",
		figOUT = lambda wildcards: wildcards.ext
	conda:
		"workflow/envs/seurat.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src[0]} {input.barcodes} {input.norm_mtx} {params.ds} {params.query} {params.ann} {params.figOUT} {output}"

### Export
rule exp_seuratobj:
	input:
		src = ["workflow/scripts/export_seurat.R", "workflow/src/10Xmat.R"],
		umi = expand("out/data2/Contrasts/{{cid}}/{{region}}/{fl}", fl=["barcodes.tsv", "features.tsv", "matrix.mtx"]),
		#NOTE: annotated barcodes replaces barcodes of norm data
		barcodes = "out/ann/Contrasts/{cid}/{region}/logNormSCT_harmony_barcodes.tsv",
		sct = expand("out/norm/Contrasts/{{cid}}/{{region}}/logNormSCT_{fl}", 
				fl=["features.tsv", "counts.mtx", "data.mtx", "scaled.tsv", "HGV.txt"]),
		harmony = expand("out/dim/Contrasts/{{cid}}/{{region}}/logNormSCT_harmony_{suffix}", suffix=["embeddings.tsv", "loadings.tsv", "projected.tsv", "stdev.txt", "nPCs.txt"]),
		graph = expand("out/clust/Contrasts/{{cid}}/{{region}}/logNormSCT_harmony_{suffix}", suffix=["GraphNN.rds", "GraphSNN.rds"])
	output:
		fl = "out/export/SeuratObject/Zhang2020_AgingArteryAtlasCynomolgusMonkey_{cid}_{region}.rds"
	params:
		GID = lambda wildcards: wildcards.cid,
		REG = lambda wildcards: wildcards.region,
		red = "harmony"

	conda:
		"workflow/envs/seurat.yaml"
	shell:
		"Rscript --vanilla "
		"{input.src[0]} {input.umi} {input.barcodes} {input.sct} {input.harmony} {input.graph} {params.GID} {params.REG} {params.red} {output.fl}"


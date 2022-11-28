import os
import pandas as pd
from argparse import ArgumentParser
from joblib import Parallel, delayed

def align(fileName, cpu, bismarkPath, output_dir, fastqPath):
    print("Step 1: Using bismark tool to align list of sequencing data to BAM file")
    with open(fileName) as fp:
        list_fileName = fp.read().split()
    def process(index, sample):
        print(f"THEARD: {index}")
        command=f"{bismarkPath}/bismark --genome {fastqPath} {sample} --output_dir {output_dir} --minimap2 --mm2_pacbio"
        print(command)
        os.system(command)
    results = Parallel(n_jobs=cpu)(delayed(process)(index, value) for index, value in enumerate(list_fileName))
    print(results)

def methyl_extract(outputDir):
    print("Step 2: Using bismark tool to extract methylation information from BAM file")
    listBam = [i for i in os.listdir(outputDir) if "bam" in i]
    print(listBam)
    def process(sample):
        command=f"{bismarkPath}/bismark_methylation_extractor --gzip --bedGraph {os.path.join(outputDir, sample)} --output_dir {outputDir}"
        print(command)
        os.system(command)

    results = Parallel(n_jobs=cpu)(delayed(process)(i) for i in listBam)
 

def bam_2_bed(outputDir):
    print("Step 3.1: Converting BAM to BED format")
    listBam = [i for i in os.listdir(outputDir) if "bam" in i]
    print(listBam)
    def process(sample):
        command=f"bedtools bamtobed -i {os.path.join(outputDir, sample)} > {os.path.join(outputDir, sample[:-4])}.bed"
        print(command)
        os.system(command)
        return f"{os.path.join(outputDir, sample[:-4])}.bed"
    results = Parallel(n_jobs=cpu)(delayed(process)(i) for i in listBam)
    return results

def mapping_CpG(outputDir):
    print("Step 3.2: Using annotatr library to map id with CpG island")
    listBed = [i for i in os.listdir(outputDir) if i.endswith(".bed")]
    def process(sample):
        command = f"Rscript --no-save script/annotation_cpg.r {os.path.join(outputDir, sample)} {os.path.join(outputDir, sample[:-4])}.csv"
        print(command)
        os.system(command)
        return f"{os.path.join(outputDir, sample[:-4])}.csv"
    results = Parallel(n_jobs=cpu)(delayed(process)(i) for i in listBed)
    return results

def combine(file_bed, file_mapping, output_dir):
    print("Step 4: Combine mapping file and methylation information file")
    offset = 100
    df_final = pd.DataFrame()
    for i in range(len(file_bed)):
        sampleName = file_bed[i].split("/")[-1][:-4]
        df_mapping = pd.read_csv(file_mapping[i])
        df_rs = pd.read_csv(os.path.join(output_dir, f'CpG_{sampleName}.txt'), skiprows=1, sep='\t', header=None)
        df_rs['unmethylated'] = (df_rs[4]=="z").astype(int)
        df_rs['methylated'] = (df_rs[4]=="Z").astype(int)
        df_mapping.set_index("name", inplace=True)
        df_rs.rename(columns={0:'name'}, inplace=True)
        tmp = df_rs.merge(df_mapping, how='inner', on='name')
        df = tmp[['annot.id', 'methylated', 'unmethylated']].groupby(by=['annot.id']).sum()
        df["betaValue"] = df['methylated']/(df['methylated'] + df['unmethylated']+offset)
        df.index = df.index.str.replace("island:", 'cpg_')
        df = df[["betaValue"]].rename(columns={"betaValue":sampleName})
        if len(df_final) == 0:
            df_final = df
        else:
            df_final = df_final.merge(df, how='inner', on='annot.id')
        df_final.to_csv(os.path.join(output_dir,"betaMatrix.csv"))
if __name__ == '__main__':
    description = 'DNAm Pipeline'
    parser = ArgumentParser(description=description, add_help=False)
    parser.add_argument('--file-name', type=str, required=True,
                        dest='fileName', help='Path file of fastq files list')
    parser.add_argument('--cpu', type=int, required=False,
                        dest='cpu', help='Parallel')
    parser.add_argument('--bismark', type=str, required=True,
                        dest='bismarkPath', help='bismarks path')
    parser.add_argument('--output-dir', type=str, required=True,
                        dest='output_dir', help='output path')
    parser.add_argument('--fastq', type=str, required=True,
                        dest='fastq', help='fastq path')

    args = parser.parse_args()
    fileName = args.fileName
    cpu = args.cpu
    bismarkPath = args.bismarkPath
    outputDir = args.output_dir
    fastqPath = args.fastq
    align(fileName, cpu, bismarkPath, outputDir, fastqPath)
    methyl_extract(outputDir)
    file_bed = bam_2_bed(outputDir)
    file_mapping = mapping_CpG(outputDir)
    combine(file_bed, file_mapping, outputDir)
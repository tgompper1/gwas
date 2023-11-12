# Tess Gompper #260947251
from scipy.stats import chi2_contingency
import allel
import pandas as pd

# 0 = reference allele, 1= a lternate allele
# disease association table looks like:
    #     | genotype 0      | genotype 1    | genotype 2
    # ----|-----------------|---------------|---------------
    # Y=0 | # with hom_ref  | # with het    | # with hom_alt
    # ----|-----------------|---------------|---------------
    # Y=1 | # with hom_ref  | # with het    | # with hom_alt
    # ----|-----------------|---------------|---------------

# Load SNPs from the VCF file
def load_vcf(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Initialize lists to store SNP genotypes and SNP IDs
    snp_genotypes = []

    for line in lines:
        if line.startswith('#'):
            continue  # Skip header lines
        fields = line.strip().split('\t')
        snp_genotypes.append(fields[6:])
    return snp_genotypes

# Load Phenotypes from the txt file
def load_phenotypes(file_path):
    arr = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
    for line in lines:
        fields = line.strip().split('\t')
        arr.append(fields[1])
    i = 0
    for p in arr:
        arr[i] = int(arr[i])
        i+=1
    return arr

if __name__ == "__main__":
    print("running...")
    # Define file paths for the VCF and phenotype data
    vcf_file = "gwas_population.vcf"
    phenotype_file = "gwas_phenotypes.txt"
    
    # Load SNP genotypes from the VCF file
    snp_genotypes = load_vcf(vcf_file)
    # Load phenotypes from the txt file
    phenotypes = load_phenotypes(phenotype_file)

    df = pd.DataFrame()
    df = pd.DataFrame(columns=['SNP_ID', 'uncorrected_pvalue', 'corrected_pvalue', 'disease_odds_het', 'disease_odds_hom', 'uncorrected_below_val'])

    ## do this for each snp
    p_value_table = []
    s = 0
    for snp in snp_genotypes:
        arr = []
        tracker=0
        for ind in snp_genotypes[s]:
            individual = ind.split('|')
            individual[0] = int(individual[0])
            individual[1] = int(individual[1])
            i = []
            i.append(individual)
            arr.append(i)
            tracker+=1
        
        genotype_array = allel.GenotypeArray(arr)

        healthy_ref, healthy_het, healthy_alt = 0, 0, 0
        diseased_ref, diseased_het, diseased_alt = 0, 0, 0
        p = 0
        for person in phenotypes:
            if person == 0:
                if genotype_array.is_hom_ref()[p][0]==True:
                    healthy_ref+=1
                elif genotype_array.is_het()[p][0]==True:
                    healthy_het+=1
                elif genotype_array.is_hom_alt()[p][0]==True:
                    healthy_alt+=1
            elif person == 1:
                if genotype_array.is_hom_ref()[p][0]==True:
                    diseased_ref+=1
                elif genotype_array.is_het()[p][0]==True:
                    diseased_het+=1
                elif genotype_array.is_hom_alt()[p][0]==True:
                    diseased_alt+=1
            p+=1
        num_skipped = 0
        disease_association_table = [[healthy_ref, healthy_het, healthy_alt],[diseased_ref, diseased_het, diseased_alt]]
        #disease_association_table_adjusted = 
        if healthy_ref == 0 or healthy_het == 0 or healthy_alt== 0  or diseased_ref==0 or diseased_het == 0 or diseased_alt== 0:
            num_skipped += 1 ## not sure about this -- the chi2 wont work with zeros so idk what to do when theres a zero 
        else:
            id = "snp"+str(s)
            res = chi2_contingency(disease_association_table)   
            p_value_u = res.pvalue
            p_value_c = res.pvalue*(10000-num_skipped)#  In the pdf he said it's product but im quite sure its division , not sure that its 3 though 
            if p_value_u < 0.05:
                below = True
            else:
                below = False
            ## Pr[Disease | individual is homo-ref]
            pr_dis_homo_ref = (diseased_ref) / (diseased_ref+healthy_ref)
            ## Calculate disease odds for heterozygous individuals
            ### Odds_ratio(het) = Pr [Disease | individual is het] / Pr [Disease | individual is homo-ref]
            odd_ratio_het = (diseased_het/(diseased_het+healthy_het)) / pr_dis_homo_ref
            ## Calculate disease odds for homozygous alternate individuals
            ### Odds_ratio(homo-alt) = Pr [Disease | individual is homo-alt] / Pr [Disease | individual is homo-ref]
            odd_ratio_homo_alt = (diseased_alt/(diseased_alt+healthy_alt)) / pr_dis_homo_ref

            #['SNP_ID', 'uncorrected_pvalue', 'corerected_pvalue', 'disease_odds_het', 'disease_odds_hom', 'uncorrected_below_val'])
            df.loc[len(df.index)] = [id, p_value_u, p_value_c, odd_ratio_het, odd_ratio_homo_alt, below] 
        s+=1  

    # count number of SNPs with an uncorrected p-value below 0.05
    count = df['uncorrected_below_val'].value_counts()[True]

    # Display SNPs with a significant associaHon after correcting for multiple 
    #   hypothesis testing using the Bonferonni correction.
    result = df.loc[df['corrected_pvalue'] <= 0.05] 
    
    print(result) 
    print("number of SNPs with an uncorrected p-value below 0.05: " + str(count))
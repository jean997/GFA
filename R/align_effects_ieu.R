
#'@title Align alleles looked up from IEU Open GWAS using ieugwasr package
#'@param phe Output of ieugwasr::phewas or ieugwasr::associations
#'@param ref Reference data frame
#'@param ref_snp_col Name of snp column in ref
#'@param ref_ea Name of effect allele column in ref
#'@param ref_nea Name of non effect allele column in ref
#'@export
align_effects_ieu <- function(phe, ref, ref_snp_col = "snp",
                          ref_ea = "A1", ref_nea = "A2"){

  ref <- ref %>%
    select(sym(ref_snp_col), sym(ref_ea), sym(ref_nea)) %>%
    rename(rsid = sym(ref_snp_col),
           A1 = sym(ref_ea),
           A2 = sym(ref_nea))
  X <- left_join(phe, ref, by="rsid")


  X <- X %>%
    mutate(ea_flip = recode(ea, "C" = "G", "T" = "A", "A" = "T", "G"="C"),
           nea_flip = recode(nea, "C" = "G", "T" = "A", "A" = "T", "G"="C"),
           sgn =  case_when(A1 == ea & A2 == nea ~ 1,
                            A1 == nea & A2 == ea ~ -1,
                            A1 == ea_flip & A2 == nea_flip ~ 1,
                            A1 == nea_flip & A2 == ea_flip ~ -1,
                            TRUE ~ 0
           ),
           beta = beta*sgn,
           eaN = case_when(sgn == 1 ~ ea,
                           sgn == -1 ~ nea,
                           TRUE ~ NA_character_),
           neaN = case_when(sgn == 1 ~ nea,
                            sgn == -1 ~ ea,
                            TRUE ~ NA_character_),
           eaf = case_when(sgn==1 ~ eaf,
                           sgn == -1 ~ 1-eaf,
                           TRUE ~ NA_real_)) %>%
    select(-A1, -A2, -ea_flip, -nea_flip, -sgn, -ea, -nea) %>%
    rename(ea = eaN, nea = neaN)
  return(X)
}

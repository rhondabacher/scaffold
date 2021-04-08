#include <Rcpp.h>
using namespace Rcpp;

CharacterVector build_ugenes(CharacterVector genes);
NumericMatrix build_count_tab(DataFrame X, CharacterVector genes);


// [[Rcpp::export]]
List get_my_tabs(List cnt_split_cell, CharacterVector genes)
{
    List L(cnt_split_cell.length());
    for (int i = 0; i < cnt_split_cell.length(); i++)
    {
        DataFrame X = cnt_split_cell[i];
        CharacterVector gene = X["Gene"];
        CharacterVector ugenes = build_ugenes(gene);
        X["ugenes"] = ugenes;
        // for (auto u : ugenes)
        //    Rcout << u << "\n";
        NumericMatrix counts_tab = build_count_tab(X, genes);
        L[i] = counts_tab;
    }

    return(L);

}

CharacterVector build_ugenes(CharacterVector genes)
{
    CharacterVector ugenes;
    for (int i = 0; i < genes.length(); i++)
    {
        String g(genes[i]);
        std::string gene_string = g.get_cstring();
        std::size_t pos = gene_string.find("@");
        std::string gene_sub = gene_string.substr(0, pos);

        String gene_sub_r(gene_sub);
        ugenes.push_back(gene_sub_r);

    }
    return(ugenes);
}

NumericMatrix build_count_tab(DataFrame X, CharacterVector genes)
{
    CharacterVector ugenes = X["ugenes"];
    CharacterVector ugenes_unique = sort_unique(ugenes);
    NumericVector all_counts = X["Count"];
    NumericVector counts_tab_tempvec(ugenes_unique.length());
    for (int i = 0; i < ugenes_unique.length(); i++)
    {
       String u = ugenes_unique[i];
       LogicalVector to_sum_logical(ugenes.length());
       for (int j = 0; j < ugenes.length(); j++)
            to_sum_logical[j] = (u == ugenes[j]);

       NumericVector to_sum = all_counts[to_sum_logical];
       float current_sum = sum(to_sum);
       counts_tab_tempvec[i] = current_sum;
    }

    counts_tab_tempvec.attr("dim") = Dimension(counts_tab_tempvec.length(), 1);
    NumericMatrix counts_tab_temp = as<NumericMatrix>(counts_tab_tempvec);
    rownames(counts_tab_temp) = ugenes_unique;


    CharacterVector zeroG = setdiff(genes, ugenes_unique);
    NumericVector zeroExprVec(zeroG.length());
    zeroExprVec.attr("dim") = Dimension(zeroExprVec.length(), 1);
    NumericMatrix zeroExpr = as<NumericMatrix>(zeroExprVec);
    rownames(zeroExpr) = zeroG;

    NumericMatrix counts_tab_temp_t = transpose(counts_tab_temp);
    NumericMatrix zeroExpr_t = transpose(zeroExpr);
    NumericMatrix counts_tab_t = cbind(counts_tab_temp_t, zeroExpr_t);
    NumericMatrix counts_tab = transpose(counts_tab_t);
    CharacterVector names = rownames(counts_tab_temp);
    for (auto n : zeroG) names.push_back(n);
    Rcout << names;

    rownames(counts_tab) = names;

    return(counts_tab);
}
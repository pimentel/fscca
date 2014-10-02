# Fast SCCA (fscca)

This is a port of NIPALS SCCA package written by Lee et. al.
\[[1](http://fafner.meb.ki.se/personal/yudpaw/?page_id=13),
[2](http://www.degruyter.com/view/j/sagmb.2011.10.issue-1/sagmb.2011.10.1.1638/sagmb.2011.10.1.1638.xml)\]
(scca). It is written for speedup in Rcpp.

Currently, only the lasso penalty function is implemented. Also, the only
cross-validation method implemented is the sequencing one-dimensional search
described in the paper.

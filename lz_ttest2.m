function t_feature = lz_ttest2(xx, yy)

xxmean = mean(xx);
yymean = mean(yy);
xxvar  = var(xx);
yyvar  = var(yy);
xxn    = length(xx);
yyn    = length(yy);
xxyy   = sqrt( ((xxn-1) * xxvar + (yyn-1) * yyvar ) / (xxn + yyn -2) );
t_feature = (xxmean - yymean) / ( xxyy * sqrt(1/xxn + 1/yyn) );
 %evaluate_ckd ecrad-rrtmg ecRad-RRTMG sw climate narrow-112 png
%evaluate_ckd RTE-RRTMGP-181204 RTE-RRTMGP lw climate narrow-256 png
%evaluate_ckd rrtmgp-nn-1.0 rrtmgp-nn-1.0 lw climate narrow-256 png
%evaluate_ckd rrtmgp-nn-1.1 NN lw climate narrow-256 png
%evaluate_ckd RTE-RRTMGP-181204 RTE-RRTMGP sw climate narrow-224 png
%evaluate_ckd rrtmgp-nn-2.0 RRTMGP-NN lw climate narrow-128 pdf
%evaluate_ckd rrtmgp-nn-2.0 RRTMGP-NN sw climate narrow-112 pdf
%evaluate_ckd rrtmgp-rr RRTMGP lw climate narrow-128 pdf
%evaluate_ckd rrtmgp-rr RRTMGP sw climate narrow-112 pdf
evaluate_ckd RTE-RRTMGP-v1.5 RRTMGP-v1.5-224 sw climate narrow-224 png
evaluate_ckd RTE-RRTMGP-v1.5 RRTMGP-v1.5-112 sw climate narrow-112 png
evaluate_ckd RTE-RRTMGP-v1.5 RRTMGP-v1.5-256 lw climate narrow-256 png
evaluate_ckd RTE-RRTMGP-v1.5 RRTMGP-v1.5-128 lw climate narrow-128 png

%metric= plot_accuracy_efficiency('RTE-RRTMGP-181204','RTE-RRTMGP','lw','climate',{'narrow'},256);
%metric= plot_accuracy_efficiency('RTE-RRTMGP-181204','RTE-RRTMGP','sw','climate',{'narrow'},224);

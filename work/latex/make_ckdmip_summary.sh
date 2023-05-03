#!/bin/bash

set -e

PLOT_DIR=/home/rd/parr/src/ckdmip-plots
DOMAIN=sw
DATASET=evaluation1
APPLICATION=climate

CKD_TOOL=ecCKD
#CKD_TOOL=ecRad-RRTMG

if [ "$DOMAIN" = lw ]
then
    XWAVE=Longwave
    xWAVE=longwave
    FLUXCODE=fluxes-4angle
else
    XWAVE=Shortwave
    xWAVE=shortwave
    FLUXCODE=fluxes
fi

if [ "$CKD_TOOL" = ecCKD ]
then
    if [ "$DOMAIN" = lw ]
    then
	BAND_STRUCTURES="fsck wide narrow"
	BAND_STRUCTURES_STR="\textit{fsck} (one full-spectrum band), \textit{wide} (5 bands) and \textit{narrow} (13 bands)"
    else
	BAND_STRUCTURES="wide narrow"
	BAND_STRUCTURES_STR="\textit{wide} (5 bands) and \textit{narrow} (13 bands)"
    fi
    APPLICATION=global-nwp
    #APPLICATION=limited-area-nwp
    #VERSION=0.5
    VERSION=0.6
    PREFIX=ecckd-${VERSION}

elif [ "$CKD_TOOL" = ecRad-RRTMG ]
then
    VERSION=1.2.0
    PREFIX=ecrad-rrtmg
    BAND_STRUCTURES="narrow"
    BAND_STRUCTURES_STR="\textit{Narrow} (13 bands)"
fi




SCENARIOS="present"

if [ "$APPLICATION" = climate ]
then
    APPSTR="Climate"
    SCENARIOS="glacialmax preindustrial present future"
elif [ "$APPLICATION" = global-nwp ]
then
    APPSTR="Global NWP"
elif [ "$APPLICATION" = limited-area-nwp ]
then
    APPSTR="Limited-Area NWP"
fi

if [ "$DATASET" = evaluation1 ]
then
    DATASETSTR="Evaluation-1"
else
    DATASETSTR=$DATASET
fi

LATEX_FILE=${PREFIX}_${DOMAIN}_${APPLICATION}_summary.tex

cat > ${LATEX_FILE} <<EOF
\documentclass{article}
\usepackage{myarticle}
\usepackage{graphicx}
\usepackage{times}
\usepackage{fancyhdr}
\usepackage{tocloft}
\usepackage{pdflscape}
\usepackage{url}
\cftsetindents{section}{2em}{5em}
\graphicspath{ {${PLOT_DIR} } }
\usepackage[multidot]{grffile}
\usepackage{rotating}
\fancyhf{}
\lhead{\bfseries\nouppercase\leftmark}
\rhead{\bfseries\thepage}
\renewcommand\thesection{Model \arabic{section}:}
\title{\Huge CKDMIP: CKD tool performance evaluation}
\pagestyle{fancy}
\begin{document}
%\setlength{cftsectionnumwidth}{6em}
\renewcommand{\cftsecfont}{\normalsize}
\setlength{\cftbeforesecskip}{0pt}
\setlength{\parindent}{0pt}
\maketitle

\begin{center}
\huge
\textmd{CKD tool:} \textbf{$CKD_TOOL version $VERSION}

\textmd{Spectral domain:} \textbf{$XWAVE}

\textmd{Application:} \textbf{$APPSTR}

\textmd{Evaluation dataset:} \textbf{$DATASETSTR}
\end{center}
\tableofcontents
\section*{Overview}
%\addcontentsline{toc}{section}{Overview}

This automatically generated document contains an evaluation of the
performance of $CKD_TOOL for generating $xWAVE correlated
\$k\$-distribution (CKD) gas-optics models targeting the application
\textit{$APPSTR}:
EOF
    if [ "$APPLICATION" = climate ]
    then
	cat >> ${LATEX_FILE} <<EOF
atmospheric heating rates are required to a minimum pressure of
0.02~hPa, and evaluation is performed for a wide range of greenhouse
gas concentrations.
EOF
    elif [ "$APPLICATION" = global-nwp ]
    then
	cat >> ${LATEX_FILE} <<EOF
atmospheric heating rates are required to a minimum pressure of
0.02~hPa, and evaluation is performed for present-day greenhouse
gas concentrations.
EOF
    else
	cat >> ${LATEX_FILE} <<EOF
atmospheric heating rates are required to a minimum pressure of
4~hPa, and evaluation is performed for present-day greenhouse
gas concentrations.
EOF
    fi

	cat >> ${LATEX_FILE} <<EOF
The evaluation dataset is \textit{$DATASETSTR} from the Correlated
K-Distribution Model Intercomparison Project
(CKDMIP)\footnote{\url{https://confluence.ecmwf.int/display/CKDMIP}}.
EOF
	if [ "$DOMAIN" = lw ]
	then
	    cat >> ${LATEX_FILE} <<EOF
Longwave radiative transfer is performed using four angles per hemisphere.
EOF
	else
	    cat >> ${LATEX_FILE} <<EOF
Shorwave radiative transfer is performed using a two-stream solver.
EOF
	fi

ACC_EFF_PLOT=${PLOT_DIR}/${PREFIX}_${DATASET}_${DOMAIN}_${APPLICATION}_accuracy-efficiency.pdf
if [ -r "$ACC_EFF_PLOT" ]
then
    cat >> ${LATEX_FILE} <<EOF
The $CKD_TOOL tool has been used to generate CKD models with the following
band structure(s): $BAND_STRUCTURES_STR. For each band structure, a number
of CKD models have been generated, characterized by the total number
of \$k\$ terms (also known as g points).
\vskip 2em

\includegraphics[width=\columnwidth]{${PLOT_DIR}/${PREFIX}_${DATASET}_${DOMAIN}_${APPLICATION}_accuracy-efficiency.pdf}
\textit{Biases and root-mean-squared errors (RMSE) in
top-of-atmosphere (TOA) upwelling irradiance and surface downwelling
irradiance, and RMSE in heating rate for two pressure ranges, for the
various band structures as a function of the total number of \$k\$
terms.
EOF
    if [ "$APPLICATION" = climate ]
    then
	if [ "$DOMAIN" = lw ]
	then
	    cat >> ${LATEX_FILE} <<EOF
It was computed from the CKDMIP scenarios 1--22, which cover climate
conditions from glacial maximum up to the worst of the CMIP6 future
scenarios, and perturbations of the individual greenhouse gases
CO\$_2\$, CH\$_4\$, N\$_2\$O, CFC11 and CFC12.}
EOF
	else
	    cat >> ${LATEX_FILE} <<EOF
It was computed from the CKDMIP scenarios 1--18, which cover climate
conditions from glacial maximum up to the worst of the CMIP6 future
scenarios, and perturbations of the individual greenhouse gases
CO\$_2\$, CH\$_4\$ and N\$_2\$O.}
EOF
	fi
    else
	cat >> ${LATEX_FILE} <<EOF
It was computed from the \`\`present-day'' CKDMIP scenario.}

EOF
    fi

fi

for BANDSTRUCT in $BAND_STRUCTURES
do
    PFILES=$(echo ${PLOT_DIR}/${PREFIX}_${DATASET}_${DOMAIN}_${APPLICATION}_${BANDSTRUCT}-*_${FLUXCODE}_present.pdf)
    NGPOINTLIST=
    for PF in $PFILES
    do
	NGPOINTLIST="$NGPOINTLIST $(basename $PF | awk -F_ '{print $5}' | awk -F- '{print $2}')"
    done
    NGPOINTLIST=$(echo $NGPOINTLIST | tr " " "\n" | sort -n)
    for NGPOINTS in $NGPOINTLIST
    do

	cat >> ${LATEX_FILE} <<EOF
\pagebreak
\section{$CKD_TOOL $APPLICATION-$BANDSTRUCT-$NGPOINTS}
%This model uses the \emph{$BANDSTRUCT} band structure with $NGPOINTS \$k\$ terms.
%\vskip 1em
%
\centerline{\includegraphics[width=0.7\columnwidth]{${PLOT_DIR}/${PREFIX}_${DOMAIN}_${APPLICATION}_${BANDSTRUCT}-${NGPOINTS}_spectral-definition.pdf}}
\textit{Illustration of the parts of the $xWAVE spectrum that
contribute to each \$k\$ term of the
$APPLICATION-$BANDSTRUCT-$NGPOINTS model.}
\vskip 1em

EOF

	if [ "$APPLICATION" = climate ]
	then
	    for SCENARIO in $SCENARIOS
	    do
		cat >> ${LATEX_FILE} <<EOF
\fbox{\includegraphics[width=0.48\columnwidth]{${PLOT_DIR}/${PREFIX}_${DATASET}_${DOMAIN}_${APPLICATION}_${BANDSTRUCT}-${NGPOINTS}_${FLUXCODE}_${SCENARIO}.pdf}}
EOF
	    done
	    if [ "$DOMAIN" = lw ]
	    then
		cat >> ${LATEX_FILE} <<EOF
\textit{Each boxed group of panels evaluate the
${APPLICATION}-${BANDSTRUCT}-${NGPOINTS} CKD model for a single CKDMIP
scenario. The left three panels in each group show the irradiances and
heating rates from the reference line-by-line calculations. The red
lines in the middle three panels show the corresponding bias in
these quantities from the CKD model. The shaded regions encompass
95\%\ of the instantaneous errors. Panels c and f
depict instantaneous errors in upwelling TOA and downwelling surface
irradiances. Error metrics are provided in the lower right.}
\vskip 2em 
EOF
	    else
		cat >> ${LATEX_FILE} <<EOF
\textit{Each boxed group of panels evaluate the
${APPLICATION}-${BANDSTRUCT}-${NGPOINTS} CKD model for a single CKDMIP
scenario. The left three panels in each group show the irradiances and
heating rates from the reference line-by-line calculations for five values of the cosine of the solar zenith angle, \$\mu_0\$. The red
lines in the middle three panels show the corresponding bias in
these quantities from the CKD model. The shaded regions encompass
95\%\ of the instantaneous errors. Panels c and f
depict instantaneous errors in upwelling TOA and downwelling surface
irradiances.
% Error metrics are provided in the lower right.
}
\vskip 2em 
EOF
	    fi

	else
	    for SCENARIO in $SCENARIOS
	    do
		cat >> ${LATEX_FILE} <<EOF

\centerline{\includegraphics[width=0.95\columnwidth]{${PLOT_DIR}/${PREFIX}_${DATASET}_${DOMAIN}_${APPLICATION}_${BANDSTRUCT}-${NGPOINTS}_${FLUXCODE}_${SCENARIO}.pdf}}
\textit{Evaluation of the ${APPLICATION}-${BANDSTRUCT}-${NGPOINTS} CKD
model for the \`\`present-day'' CKDMIP scenario. The left three panels
show the irradiances and heating rates from the reference line-by-line
calculations. The red lines in the middle three panels show the
corresponding bias in these quantities from the CKD model. The shaded
regions encompass 95\%\ of the instantaneous errors. Panels c and f
depict instantaneous errors in upwelling TOA and downwelling surface
irradiances. Error metrics are provided in the lower right.}
\vskip
2em


EOF
	    done
	fi

	if [ "$BANDSTRUCT" = wide ]
	then
	    cat >> ${LATEX_FILE} <<EOF
\includegraphics[width=\textwidth]{${PLOT_DIR}/${PREFIX}_${DATASET}_${DOMAIN}_${APPLICATION}_${BANDSTRUCT}-${NGPOINTS}_bands_present.pdf}
\textit{Evaluation of irradiances and heating rates for the broadband
(leftmost column of panels) and the 5 wide $xWAVE bands (other
panels) of the ${APPLICATION}-${BANDSTRUCT}-$NGPOINTS CKD model. The
black dashed and red solid lines correspond to the average of the 50
profiles for the \`\`present-day'' scenario, while the shaded regions
encompass 95\%\ of the error.}  \vskip 2em

EOF
	elif [ "$BANDSTRUCT" = narrow ]
	then
	    cat >> ${LATEX_FILE} <<EOF
\begin{landscape}
%\rotatebox{90}{
\includegraphics[width=24cm]{${PLOT_DIR}/${PREFIX}_${DATASET}_${DOMAIN}_${APPLICATION}_${BANDSTRUCT}-${NGPOINTS}_bands_present.pdf}
%}

\textit{Evaluation of irradiances and heating rates for the broadband
(leftmost column of panels) and the 13 narrow $xWAVE bands (other panels) of
the ${APPLICATION}-${BANDSTRUCT}-$NGPOINTS CKD model. The black
dashed and red solid lines correspond to the average of the 50
profiles for the \`\`present-day'' scenario, while the shaded regions
encompass 95\%\ of the error.}
\vskip 2em
\end{landscape}
EOF

	fi
	if [ "$APPLICATION" = climate ]
	then
	    cat >> ${LATEX_FILE} <<EOF
\centerline{\includegraphics[width=\columnwidth]{${PLOT_DIR}/${PREFIX}_${DATASET}_${DOMAIN}_${APPLICATION}_${BANDSTRUCT}-${NGPOINTS}_forcing.pdf}}
\textit{Comparison of reference line-by-line and calculations by the
${APPLICATION}-${BANDSTRUCT}-${NGPOINTS} model of the instantaneous
clear-sky radiative forcing from perturbing each of the five
well-mixed greenhouse gases from their present-day values, at (top
row) top-of-atmosphere and (middle row) surface, averaged over the 50
profiles of the $DATASETSTR dataset. The bottom row shows the mean change
to heating rate resulting from perturbing the concentration of a gas
from its present-day value to either the maximum or minimum value in
the range for that gas.}
\vskip 2em
EOF

	    if [ "$DOMAIN" = lw ]
	    then
		cat >> ${LATEX_FILE} <<EOF
\centerline{\includegraphics[width=0.75\columnwidth]{${PLOT_DIR}/${PREFIX}_${DATASET}_${DOMAIN}_${APPLICATION}_${BANDSTRUCT}-${NGPOINTS}_overlap.pdf}}
\textit{Evaluation of the representation of spectral overlap of
CO\$_2\$, CH\$_4$ and N\$_2\$O by the
${APPLICATION}-${BANDSTRUCT}-${NGPOINTS} CKD model. In each panel, the
abscissa shows the TOA radiative forcing from perturbing a gas to
either its climatic minimum or maximum value. These radiative forcings
are computed keeping the concentrations of all other well-mixed gases
at their present-day values, except for the gas on the ordinate, which is
perturbed to its own climatic minimum or maximum values.}
EOF
	    fi
	fi

    done

done

	cat >> ${LATEX_FILE} <<EOF
\end{document}
EOF

pdflatex ${LATEX_FILE}
pdflatex ${LATEX_FILE}

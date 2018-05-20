(* ::Package:: *)

(* Mathematica Package *)
(* :Context: ChainStat` *)
(* :Author: Yi Wang *)
(* :Date: 6/1/14 *)
(* The package is released under GPLv3. Report bug to tririverwangyi@gmail.com *)

BeginPackage["ChainStat`"]
S2P::usage = "S2P[dervation] gives probability corresponding to derivation. For example, S2P[1.0] = 0.682689";
P2S::usage = "P2S[probability] gives standard derivation. For example, P2S[0.682689] = 1";
CntList::usage = "CntList[data, \"Operation\"->Opt, \"NBin\"->Num] bins data of any dimensions.
	Here data = {{x1, y1, ...}, {x2, y2, ...}, ...}.
	The return value is {{xbin1, ybin1, ..., count1}, {xbin2, ybin2, ..., count2}, ...},
	where {{xbin1, ybin1, ...}, {xbin2, ybin2, ...}} is a flattened grid of size nBin in each dimensions;
	and {count1, count2, ...} is number count, or Opt acting on the number count if \"Operation\" specified.";
HistList::usage = "HistList[data, \"NBin\"->Num] is CntList, with \"Operation\" as to calculate the CDF.";

$Chain::usage = "Chain in operation."
SmoothContourPlot::usage = "SmoothContourPlot[data] is similar to ListContourPlot, while smooths the result.";
SmoothDensityPlot::usage = "SmoothDensityPlot[data] is similar to ListDensityPlot, while smooths the result."
LoadChain::usage = "LoadChain[filename, burnInRate] reads MCMC chains, which returns $Chain.
	$Chain contains {names, texNames, chain}, where
	chain contains {numberOfStay, chi2, value for vars in names}.
	A list of filenames can be given and in this case those chains are combined (assuming the same parameters).";
CondChain::usage = "CondChain[condition] produces new chainInfo from $Chain, subject to condition.";
AddToChain::usage = "AddToChain[varName, texName, expression] add a column to chain.
  Example: AddToChain[\"nsMinus1\", \"n_s-1\", \"ns\"-1] adds nsMinus1 to chain.";
Sel::usage = "Sel[{var1, var2, ...}, condition] selects samples of variables, subject to optional condition.";
BestFits::usage = "BestFits[n] prints the first n best fit points."
Cnt::usage = "Cot[{var1, var2, ...}, cond] counts number of samples on binned grads, subject to optional condition.";
PlotDensity::usage = "PlotDensity[{var1, var2}] performs density plot of the two variables' number counting.";
PlotContour::usage = "PlotContour[{var1, var2}] performs contour plot of the two variables' CDF.";
PlotSample::usage = "PlotSample[{var1, var2}] performs direct plot of the two variables' data points.";
Hist::usage = "Cot[{var1, var2, ...}, cond] calculates CDF on binned grads, subject to optional condition.";
Stat::usage = "Stat[var] or Stat[{var1, var2}] do 1d or 2d statistics and plots.";
HistNBins::usage = "Number of bins in 1d statistics.";
PtPlot::usage = "PtPlot[{var1, var2}, cond] plots points on two dimensions.";
TrianglePlot::usage = "TrianglePlot[{var1, var2, var3, ...}] performs triangle plots.";

Begin["`Private`"]

(* ****************************************************************************************************************** *)
(* Utilities *)

P2SmaxS=10.0; (* to prevent zero probability and infinity sigma from binning *)
P2S = Sqrt[2.]InverseErf[#] /. \[Infinity]->p2SmaxS &;
S2P = Erf[#/Sqrt[2.]]&;

(* ****************************************************************************************************************** *)
(* Binning of list *)

Options[CntListM] = Options[CntListW] = {"Operation"->Identity, "NBin"->30};
CntList=CntListW;

(* M version uses BinCounts. However, BinCounts needs to be applied to grids which is slow and miss points. Otherwise boundary is not regular. *)
CntListM[data_, OptionsPattern[]]:= Module[
  {min, max, binCoord, binCount, nBin=OptionValue["NBin"], operation=OptionValue["Operation"]},
  {min, max} = Transpose[{Min@#, Max@#}& /@ Transpose[data]];
  binCoord = # * (max-min) / nBin + min & /@ Tuples[Range[nBin]-0.5, Length@data[[1]]];
  binCount = Flatten @ BinCounts[data, Apply[Sequence,{Range[#1, #2, (#2-#1)/nBin]}& @@@ Transpose[{min,max}]] ];
  MapThread[Append[#1,#2]&,{binCoord, operation@binCount}] ]

(* Homemade BinCounts. Faster and does not miss points. *)
CntListW[data_, OptionsPattern[]]:=Module[{tmin, tmax, nBin=OptionValue["NBin"], operation=OptionValue["Operation"], sumCnt, sumOne, f, res},
  {tmax, tmin} = Transpose[{Max@#, Min@#} & /@ Transpose[data]];
  sumCnt = Total[Apply[f, Ceiling[$MachineEpsilon + nBin (# - tmin)/(tmax - tmin)]] & /@ data];
  sumOne = Total[f @@@ Tuples[Range@nBin, Length@tmax]];
  res = Append[({##}[[;;-2]]-0.5) (tmax-tmin)/nBin + tmin, {##}[[-1]]] & @@@ (Apply[List, sumCnt+sumOne] /. n_. f[a__] :> {a,n-1});
  Transpose @ Append[Transpose[res][[;;-2]], operation@Transpose[res][[-1]]] ]

CalcCDF[binCount_]:= Module[{sortFlat, cdf},
  sortFlat = Sort @ binCount;
  cdf = Accumulate[sortFlat] / N@Total[sortFlat];
  binCount /. Dispatch[MapThread[Rule,{sortFlat, cdf}]] (* this is CDF *)]

HistList[data_, opt:OptionsPattern[CntList]]:= CntList[data, "Operation"->CalcCDF, opt];

(* ****************************************************************************************************************** *)
(* Plotting *)

hasLMRfont = !FreeQ[FE`Evaluate@FEPrivate`GetPopupList@"MenuListFonts","Latin Modern Roman"];
PlotOptions = {Frame->True, PlotRange->All, Axes->False, ImageSize->Large, AspectRatio->1,
  ImagePadding->{{150,30},{150,30}},
  BaseStyle->{FontFamily->If[hasLMRfont, "Latin Modern Roman", "Times"], FontSize->24}};

(* Loess and LoessFit are based on code from Rahul Narain *)
Loess[nearest_, k_][x_, y_] := Module[{nearPts, d, u, v},
  nearPts = nearest[{x, y}, k];
  d = EuclideanDistance[{x, y}, Most[#]] & /@ nearPts;
  d = d/Max[d];
  LinearModelFit[nearPts, {u, v, u^2, v^2, u v}, {u, v}, Weights -> (1 - d^3)^3][x,y]];

Options[LoessFit] = {"Smoothing"->2};
LoessFit[data_List, OptionsPattern[]]:= Module[{sd, scaledData, nearest, nNear},
  sd = Most@StandardDeviation@data;
  scaledData = Transpose @ Append[Most@Transpose@data/sd, Last@Transpose@data];
  nearest = Nearest[scaledData /. {x_, y_, z_} :> ({x, y} -> {x, y, z})];
  nNear = Floor[OptionValue["Smoothing"]*Sqrt@Length@data];
  Loess[nearest, nNear][#1/sd[[1]], #2/sd[[2]]]& ]

SmoothContourPlot[data_List, opts:OptionsPattern[ContourPlot]]:= Module[{fit, xmin, xmax, ymin, ymax},
  fit = LoessFit[data];
  {{xmin, xmax}, {ymin, ymax}} = {Min[#], Max[#]} & /@ Most[Transpose[data]];
  ContourPlot[fit[x, y], {x, xmin, xmax}, {y, ymin, ymax}, opts, Contours->1-{S2P[1], S2P[2]}, ContourShading->None,
    ContourStyle->{Directive[Darker@Blue,Thick], Directive[Lighter@Blue,Thick, Dashed]}, Evaluate[Sequence@@PlotOptions]]];

SmoothDensityPlot[data_List, opts:OptionsPattern[DensityPlot]]:= Module[{fit, xmin, xmax, ymin, ymax, colorFunc},
  fit = LoessFit[data];
  {{xmin, xmax}, {ymin, ymax}} = {Min[#], Max[#]} & /@ Most[Transpose[data]];
  colorFunc = Function[z, RGBColor[1-z, 1-z, 1-z]];
  DensityPlot[fit[x, y], {x, xmin, xmax}, {y, ymin, ymax}, opts, PlotRange->All,
    ColorFunction->colorFunc, Evaluate[Sequence@@PlotOptions]]]

(* PlotDensity:= PlotDensityM; *)(* This is the faster and the smoother between the two, but has less control. *)
PlotDensity:= PlotDensityM;
PlotDensityW[{var1_String, var2_String}, cond_:True, opt:OptionsPattern[DensityPlot]] /;IsntRule[cond] :=
  SmoothDensityPlot[Cnt[{var1, var2}, cond], opt, PlotRange->All, FrameLabel->{StyledName@var1, StyledName@var2}];

PlotDensityM[{var1_String, var2_String}, cond_:True, opt:OptionsPattern[SmoothDensityHistogram]] /;IsntRule[cond] :=
  SmoothDensityHistogram[Sel[{var1, var2}, cond], opt, ColorFunction->(GrayLevel[1-#]&), PlotRange->Full,
    FrameLabel->{StyledName@var1, StyledName@var2}, Evaluate[Sequence@@PlotOptions]];

PlotContour[{var1_String, var2_String}, cond_:True, opt:OptionsPattern[ListContourPlot]] /;IsntRule[cond] :=
  SmoothContourPlot[Hist[{var1, var2}, cond], opt, FrameLabel->{StyledName@var1, StyledName@var2}];

PlotSample[{var1_String, var2_String}, cond_:True, opt:OptionsPattern[ListPlot]] /;IsntRule[cond] := ListPlot[Sel[{var1, var2}, cond],
  opt, PlotStyle -> Directive[Black, Opacity@0.05], PlotRange->All,
  FrameLabel->{StyledName@var1, StyledName@var2}, PlotOptions];
                
TrianglePlot[vars_List]:=Grid@Table[If[v1=!=v2,If[!OrderedQ[{v1,v2}],PlotDensity[{v2,v1}]],Stat[v1]],{v1,vars},{v2,vars}];    

(* ****************************************************************************************************************** *)
(* Load chain *)

Options[LoadChain] = {"BurnIn"->0.3};
LoadChain[fn_String, OptionsPattern[]]:= Module[{nameFile, names, chainFiles, chainList, burnIn=OptionValue["BurnIn"]},
  nameFile = First @ FileNames[fn<>".paramnames"];
  names = Map[StringTrim, Import[nameFile, "TSV"], -1];
  chainFiles = FileNames[(fn<>"_") ~~DigitCharacter..~~ ".txt"];
  chainList = Function[f, #[[Floor[burnIn*Length@#];;]]& @ ReadList[f,Table[Number,{Length@names+2}]]] /@ chainFiles;
  $Chain = Append[Transpose@names, Join@@chainList]]

LoadChain::difvar = "Warning: Chains have different variables. Combination may be wrong!";
LoadChain[flist_List, opt:OptionsPattern[]]:= Module[{infos = LoadChain[#, opt]& /@ flist},
  If[Length@DeleteDuplicates[Most/@infos] =!= 1, Message@LoadChain::difvar];
  $Chain = Append[infos[[1,;;-2]], Join@@Last/@infos]];

(* ****************************************************************************************************************** *)
(* Parameter name and position *)

FindName::nofnd = "Error: Name `1` not found in `2`.";
FindName::nouqe = "Error: Name `1` not unique in `2`.";
FindName[name_String]:= Module[{pos},
  pos = Position[$Chain[[1]],name];
  If[Length@pos==0, Message[FindName::nofnd, name, $Chain[[1]]]; Return[]];
  If[Length@pos>1, Message[FindName::nouqe, name, $Chain[[1]]]; Return[]];
  pos[[1,1]]]

ChainPos[name_String]:= If[name === "likelihood", 2, FindName[name]+2];

TexName["likelihood"]:="-log(likelihood)";
TexName[name_String]:= ToExpression[$Chain[[2, FindName[name]]], TeXForm, HoldForm];

StyledName[name_String]:= Style[TexName[name], 36]

(* ****************************************************************************************************************** *)
(* Selection from chain *)

ExprOnChain::usage = "ExprOnChain[expr] converts strings in expr into pure function evolving chain positions.
  For example, ExprOnChain[\"likelihood\"<4900] returns a pure function #[[2]] < 4900 &";
ExprOnChain[expr_]:= Module[{vars = Cases[{expr}, _String, Infinity], slt},
  (slottedExpr&) /. slottedExpr -> expr /. MapThread[Rule, {vars, slt /@ ChainPos /@ vars}] /.  slt -> Function[x, #[[x]]] ] ;

CondChain[cond_]:= Append[$Chain[[;;-2]], Select[$Chain[[-1]], ExprOnChain@cond]];

AddToChain[vName_, texName_, exprRaw_]:= Module[{expr, newColumn},
  expr = ExprOnChain @ exprRaw;
  newColumn = expr /@ $Chain[[3]];
  $Chain = {Append[$Chain[[1]], vName], Append[$Chain[[2]], texName], Transpose@Append[Transpose[$Chain[[3]]], newColumn]}];

Sel[vars_List]:= Flatten[#, 1]& @ (ConstantArray[{##2}, Round@#1]& (* expand numberOfStay *) @@@
    $Chain[[3]][[All, Prepend[ChainPos /@ vars, 1] (* also get numberOfStay *) ]]);
Sel[var_String]:= Flatten @ Sel @ {var};
Sel[var_, cond_]:= Block[{$Chain = CondChain[cond]}, Sel[var]];

BestFits[nPts_:5]:= Module[{names, fits},
  names = {Prepend[$Chain[[1]], "likelihood"], Prepend[TexName /@ $Chain[[1]], TexName["likelihood"]]};
  fits = SortBy[$Chain[[3]], #[[2]] &][[;; nPts, 2 ;;]];
  MatrixForm @ Join[names, fits]]

(* ****************************************************************************************************************** *)
(* Chain statistics *)

IsntRule[f_] := Head[f]=!=List && Head[f]=!=Rule;

Cnt[vars_List, cond_:True, opt:OptionsPattern[CntList]] /;IsntRule[cond] := CntList[Sel[vars, cond], opt];
Hist[vars_List, cond_:True, opt:OptionsPattern[CntList]] /;IsntRule[cond] := Cnt[vars, cond, "Operation"->CalcCDF, opt]

HistNBins = 20;
Stat[var_String, cond_:True, opt:OptionsPattern[Histogram]] /;IsntRule[cond] := Module[{selected = Sel[{var}, cond]},
  Print[TexName@var, " = ", Mean@Flatten@selected, " \[PlusMinus] ", StandardDeviation@Flatten@selected, " (mean, std var)"];
  (*Print["\[Sigma] contours: ", {#1,p2\[Sigma][1-#2]}& @@@ HistList[selected,30]];*)
  Histogram[Flatten@selected, HistNBins, "PDF", opt, FrameLabel->{StyledName@var, None},
    FrameTicks->{{None,None},{Automatic,Automatic}}, ChartStyle->Directive[EdgeForm[None], Lighter@Gray], PlotOptions]];

Stat[{var1_String, var2_String}, etc___] := PlotContour[{var1, var2}, etc];

(* ****************************************************************************************************************** *)
(* Setup of environment *)

SetDirectory[If[$Notebooks && Quiet[Check[NotebookDirectory[],False]]=!=False (* If exists saved nb *),
  NotebookDirectory[], Import["!pwd","Text"]]];

Print["ChainStat, by Yi Wang (2014, 2018, GPLv3). Report bug to tririverwangyi@gmail.com\nPlease check with independent methods before using the results."];

End[]
EndPackage[]

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
LoadChain::usage = "LoadChain[filename, burnInRate] reads MCMC chains, which returns $Chain.
	$Chain contains {names, texNames, chain}, where
	chain contains {numberOfStay, chi2, value for vars in names}.
	A list of filenames can be given and in this case those chains are combined (assuming the same parameters).";
CondChain::usage = "CondChain[condition] produces new chainInfo from $Chain, subject to condition.";
Sel::usage = "Sel[{var1, var2, ...}, condition] selects samples of variables, subject to optional condition.";
Cnt::usage = "Cot[{var1, var2, ...}, cond] counts number of samples on binned grads, subject to optional condition.";
Hist::usage = "Cot[{var1, var2, ...}, cond] calculates CDF on binned grads, subject to optional condition.";
Stat::usage = "Stat[var] or Stat[{var1, var2}] do 1d or 2d statistics and plots.";
PtPlot::usage = "PtPlot[{var1, var2}, cond] plots points on two dimensions.";

Begin["`Private`"]

(* ****************************************************************************************************************** *)
(* Utilities *)

P2SmaxS=10.0; (* to prevent zero probability and infinity sigma from binning *)
P2S = Sqrt[2.]InverseErf[#] /. \[Infinity]->p2SmaxS &;
S2P = Erf[#/Sqrt[2.]]&;

(* ****************************************************************************************************************** *)
(* Binning of list *)

Options[CntList] = {"Operation"->Identity, "NBin"->30};
CntList[data_, OptionsPattern[]]:= Module[
  {min, max, binCoord, binCount, nBin=OptionValue["NBin"], operation=OptionValue["Operation"]},
  {min, max} = Transpose[{Min@#, Max@#}& /@ Transpose[data]];
  binCoord = # * (max-min) / nBin + min & /@ Tuples[Range[nBin]-0.5, Length@data[[1]]];
  binCount = Flatten @ BinCounts[data, Apply[Sequence,{Range[#1, #2, (#2-#1)/nBin]}& @@@ Transpose[{min,max}]] ];
  MapThread[Append[#1,#2]&,{binCoord, operation@binCount}] ]

CalcCDF[binCount_]:= Module[{sortFlat, cdf},
  sortFlat = Sort @ binCount;
  cdf = Accumulate[sortFlat] / N@Total[sortFlat];
  binCount /. Dispatch[MapThread[Rule,{sortFlat, cdf}]] (* this is CDF *)]

HistList[data_, opt:OptionsPattern[CntList]]:= CntList[data, "Operation"->CalcCDF, opt];

(* ****************************************************************************************************************** *)
(* Plotting *)

hasLMRfont = !FreeQ[FE`Evaluate@FEPrivate`GetPopupList@"MenuListFonts","Latin Modern Roman"];
PlotOptions = {Frame->True, PlotRange->All, Axes->False, ImageSize->Large, AspectRatio->1,
  ImagePadding->{{150,20},{150,20}},
  BaseStyle->{FontFamily->If[hasLMRfont, "Latin Modern Roman", "Times"], FontSize->24}};

(* Code of loess and smoothContourPlot mainly comes from Rahul Narain *)
Loess[nearest_, k_][x_, y_] := Module[{nearPts, d, u, v},
  nearPts = nearest[{x, y}, k];
  d = EuclideanDistance[{x, y}, Most[#]] & /@ nearPts;
  d = d/Max[d];
  LinearModelFit[nearPts, {u, v, u^2, v^2, u v}, {u, v}, Weights -> (1 - d^3)^3][x,y]];

Options[SmoothContourPlot] = Join[{"Smoothing"->2, Contours->1-{S2P[1], S2P[2]}, ContourShading->None,
  ContourStyle->{Directive[Darker@Blue,Thick], Directive[Lighter@Blue,Thick, Dashed]}}, PlotOptions];
SmoothContourPlot[data_List, opts:OptionsPattern[{SmoothContourPlot, ListContourPlot}]]:= Module[
  {sd, scalex, scaley, scaledData, nearest, nNear, fit, xmin, xmax, ymin, ymax,
    pltOpts = Sequence @@ FilterRules[{opts}~Join~Options[SmoothContourPlot], Options[ContourPlot]] },
  sd = StandardDeviation[data];
  {scalex, scaley} = 1/Most[sd];
  scaledData = {scalex, scaley, 1} # & /@ data;
  nearest = Nearest[scaledData /. {x_, y_, z_} :> ({x, y} -> {x, y, z})];
  nNear = Floor[OptionValue["Smoothing"]*Sqrt@Length@data];
  fit = Loess[nearest, nNear];
  {{xmin, xmax}, {ymin, ymax}} = {Min[#], Max[#]} & /@ Most[Transpose[data]];
  ContourPlot[fit[scalex x, scaley y], {x, xmin, xmax}, {y, ymin, ymax}, Evaluate@pltOpts] ]

(* ****************************************************************************************************************** *)
(* Load chain *)

Options[LoadChain] = {"BurnIn"->0.3};
LoadChain[fn_String, OptionsPattern[]]:= Module[{names, chainFiles, chain, totChain={}, burnIn=OptionValue["BurnIn"]},
  names = Import[fn<>".paramnames"];
  (* for simplicity, assuming #chains \[LessEqual] 9. Otherwise, only first 9 are read. *)
  chainFiles = Import["!ls "<>fn<>"_?.txt","List"] ;
  Do[
    chain = ReadList[f,Table[Number,{Length@names+2}]];
    totChain = Join[totChain,chain[[Floor[burnIn*Length@chain];;]]]
    ,{f,chainFiles}];
  $Chain = Append[names\[Transpose],totChain]];

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

CondChain[condRaw_]:= Module[{vars, cond},
  vars = Cases[condRaw, _String, Infinity];
  cond[point_]:= condRaw /. MapThread[Rule, {vars, point[[ChainPos@#]]& /@ vars}];
  Append[$Chain[[;;-2]], Select[$Chain[[-1]], cond]]];

Sel[vars_List]:= Flatten[#, 1]& @ (ConstantArray[{##2}, Round@#1]& (* expand numberOfStay *) @@@
    $Chain[[3]][[All, Prepend[ChainPos /@ vars, 1] (* also get numberOfStay *) ]]);
Sel[var_String]:= Flatten @ Sel @ {var};
Sel[var_, cond_]:= Block[{$Chain = CondChain[cond]}, Sel[var]];

(* ****************************************************************************************************************** *)
(* Chain statistics *)

IsntRule[f_] := Head[f]=!=List && Head[f]=!=Rule;

Cnt[vars_List, cond_:True, opt:OptionsPattern[CntList]] /;IsntRule[cond] := CntList[Sel[vars, cond], opt];
Hist[vars_List, cond_:True, opt:OptionsPattern[CntList]] /;IsntRule[cond] := Cnt[vars, cond, "Operation"->CalcCDF, opt]

Stat[var_String, cond_:True, opt:OptionsPattern[Histogram]] /;IsntRule[cond] := Module[{selected = Sel[{var}, cond]},
  Print[TexName@var, " = ", Mean@Flatten@selected, " \[PlusMinus] ", StandardDeviation@Flatten@selected, " (mean, std var)"];
  (*Print["\[Sigma] contours: ", {#1,p2\[Sigma][1-#2]}& @@@ HistList[selected,30]];*)
  Histogram[Flatten@selected, FrameLabel->{StyledName@var, None}, FrameTicks->{{None,None},{Automatic,Automatic}},
    ChartStyle->Directive[EdgeForm[None], Lighter@Gray], opt, PlotOptions]];

Stat[{var1_String, var2_String}, cond_:True, opt:OptionsPattern[{SmoothContourPlot, ListContourPlot}]] /;IsntRule[cond] :=
    SmoothContourPlot[Hist[{var1, var2}, cond], FrameLabel->{StyledName@var1, Column[{StyledName@var2,""}]}, opt ]

PtPlot[{var1_String, var2_String}, cond_:True] /;IsntRule[cond] := ListPlot[Sel[{var1, var2}, cond],
  PlotStyle -> Directive[Black, Opacity@0.05], FrameLabel->{StyledName@var1, StyledName@var2}, PlotOptions];

(* ****************************************************************************************************************** *)
(* Setup of environment *)

SetDirectory[If[$Notebooks && Quiet[Check[NotebookDirectory[],False]]=!=False (* If exists saved nb *),
  NotebookDirectory[], Import["!pwd","Text"]]];

End[]
EndPackage[]
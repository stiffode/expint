#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input::Initialization:: *)
(*This package turns Mathematica input into Tex code*)Needs["CellsToTeX`"]
(*NCAlgebra to allow for non commutative expr*)
(*Import["https://raw.githubusercontent.com/NCAlgebra/NC/master/NCExtras/NCWebInstall.m",CharacterEncoding->"UTF8"];*)
(*Needs["NCAlgebra`"]*)
<<NC`;
<<NCAlgebra`;


(* ::Input::Initialization:: *)
(* To clear all variables *)
Clear["Global`*"]

(*Modify the number of stages*)
stages=4;
s=stages;
(*Building the Butcher Tableau and vector of c_i*)
a=Table[Symbol["a"<>ToString[i]<>ToString[j]],{i,0,stages-1},{j,0,i}]
c=Table[Symbol["c"<>ToString[i]],{i,0,stages}]
(*Evaluating each phi function at each c_i*)
phi1c=Table[phi1[-c[[i]]*L*h],{i,1,stages}]
phi2c=Table[phi2[-c[[i]]*L*h],{i,1,stages}]
phi3c=Table[phi3[-c[[i]]*L*h],{i,1,stages}]

(* Set the NonCommutative variables *)
(* We only need to make L non-commutative and
every function of L is automatically non-commutative (i.e. Phi functions) *)
SetNonCommutative[J,K,L,a];


(* ::Input::Initialization:: *)
(* Specifying c_i and evaluating phi functions at 0 *)
c0=0;c1=1/2;c2=3/4;
c3=1;
phi1[0]=1;phi2[0]=1/2;phi3[0]=1/(3!);phi4[0]=1/(4!);phi5[0]=1/(5!);


(* ::Input::Initialization:: *)
(* Joining all the phi functions as lists based on whether they are phi1, phi2 etc *)
phiFunctions = {Join[phi1c,{phi1[-L*h]}],Join[phi2c,{phi2[-L*h]}],Join[phi3c,{phi3[-L*h]}]};
(* VARIABLES TO SOLVE FOR *)
vars = {};
(* When some c_i are the same *)
optimizer = {};
(* Looping over weights a_ij using i and j *)
For[i=0,i<stages,i++, 
	tempListA = {};
	For[j=0,j<i+1,j++,
		tempA = 0;
		(* Looping over phi_m evaluated at c_n using m and n *)
		For[m=1,m<Length[phiFunctions]+1,m++,
		optimizer = {};
	(* n starts at 2 to avoid having phi_m at c_0 since c_0=0 *)			For[n=2,n<Length[phiFunctions[[m]]]+1,n++,
				(*Forming the linear combinations of phi for a_ij*)
If[MemberQ[optimizer,c[[n]]],,	tempA+=Symbol["a"<>ToString[i]<>ToString[j]<>ToString[m]<>ToString[n-1]]*phiFunctions[[m,n]];
SetCommutative[Symbol["a"<>ToString[i]<>ToString[j]<>ToString[m]<>ToString[n-1]]];
AppendTo[vars,Symbol["a"<>ToString[i]<>ToString[j]<>ToString[m]<>ToString[n-1]]];
AppendTo[optimizer,c[[n]]];
];
			];
	  ];
		AppendTo[tempListA,tempA];
	];
	(*Assigning the linear combinations to weights a_ij*)
	(*Crucial step so we can use them as variables later*)
	MapThread[Set,{a[[i+1]],tempListA}];
]


(* ::Input::Initialization:: *)
(* To get all coefficients in front of any phi_i eval at each c_j *)
(* We do not include phi_i at c_0 since c_0=0 *)
phiFunctionsAllTogetherNoc0 = {Join[Delete[phi1c,1],{phi1[-L*h]}],Join[Delete[phi2c,1],{phi2[-L*h]}],Join[Delete[phi3c,1],{phi3[-L*h]}]};
phiFlattened = Union[Flatten[phiFunctionsAllTogetherNoc0]];
(*To gather all possible terms involving the operators J and K *)
terms = {};
termsJ = {J};
termsK={K};
termsJJ={J**J};
For[i = 1, i<Length[phiFlattened]+1,i++,
		AppendTo[termsJ,J**phiFlattened[[i]]];
AppendTo[termsK,K**phiFlattened[[i]]];
	For[j=1,j<Length[phiFlattened]+1,j++,
		AppendTo[termsJ,phiFlattened[[i]]**J**phiFlattened[[j]]];
AppendTo[termsK,phiFlattened[[i]]**K**phiFlattened[[j]]];
	];
];
For[i = 1, i<Length[termsJ]+1,i++,
	For[j=1,j<Length[phiFlattened]+1,j++,
		AppendTo[termsJJ,termsJ[[i]]**J**phiFlattened[[j]]];
];
];

termTypes = Join[phiFlattened,termsJ,termsJJ,termsK];
(* To form a list with all Coefs for all types of terms *)
allCoefs = {};


(* ::Input::Initialization:: *)
(* To append coefs from order1lower *)
(* Equal->Subtract gets all terms on one side*)
order1lower = {Sum[a[[s]][[i]],{i,stages}]==phi1[-L*h]};
allorder1lower = order1lower;
allorder1Coeflower = CoefficientRules[allorder1lower/.Equal->Subtract,termTypes];
For[i = 1, i<Length[allorder1Coeflower]+1,i++,
	For[j=1,j<Length[allorder1Coeflower[[i]]]+1,j++,
		AppendTo[allCoefs,allorder1Coeflower[[i,j,2]]];
	];
];
Print[Length[allCoefs]]


(* ::Input::Initialization:: *)
(* Saving a_(s-1)i before evaluating their phi functions at 0 *)
asoriginal = a[[s]];
a30=a30/.L->0;
a31=a31/.L->0;
a32=a32/.L->0;
a33=a33/.L->0;
(* To append coefs from order2lower *)
(*   Applying conditions for order 3 holding weakly:            *) 
(*           a_(s-2)i(0) instead of a_(s-2)i(-hL)                         *)
(*  Always reread Theorem 4.7 in the Ho&Os05 paper when applying this step, "holding weakly" has a very particular meaning *)
order2alower = {Sum[a[[s]][[i]]*c[[i]],{i,stages}] == phi2[-L*h]};
order2blower =  Table[Sum[a[[i]][[k]],{k,1,i}]==phi1c[[i+1]]*c[[i+1]],{i,1,stages-1}];
(* The following line takes care of psi_p(0)=0 in "holding weakly" *)
order2alower = order2alower/.{L->0};
allorder2lower = Join[order2alower,order2blower];
allorder2Coeflower = CoefficientRules[allorder2lower/.Equal->Subtract,termTypes];
For[i = 1, i<Length[allorder2Coeflower]+1,i++,
	For[j=1,j<Length[allorder2Coeflower[[i]]]+1,j++,
		AppendTo[allCoefs,allorder2Coeflower[[i,j,2]]];
	];
];
Print[Length[allCoefs]]


(* ::Input::Initialization:: *)
(* To append coefs from order1higher *)
a[[s]] = asoriginal;
order1higher = {Sum[a[[s-1]][[i]],{i,stages-1}]==phi1[-L*h]};
allorder1higher = order1higher;
allorder1higherCoef = CoefficientRules[allorder1higher/.Equal->Subtract,termTypes];
For[i = 1, i<Length[allorder1higherCoef]+1,i++,
	For[j=1,j<Length[allorder1higherCoef[[i]]]+1,j++,
		AppendTo[allCoefs,allorder1higherCoef[[i,j,2]]];
	];
];
Length[allCoefs]


(* ::Input::Initialization:: *)
(* To append coefs from order2higher *)
order2ahigher = {Sum[a[[s-1]][[i]]*c[[i]],{i,stages-1}] == phi2[-L*h]};
order2bhigher = Table[Sum[a[[i]][[k]],{k,1,i}]==phi1c[[i+1]]*c[[i+1]],{i,1,stages-2}];
allorder2higher = Join[order2ahigher,order2bhigher];
allorder2higherCoef = CoefficientRules[allorder2higher/.Equal->Subtract,termTypes];
For[i = 1, i<Length[allorder2higherCoef]+1,i++,
	For[j=1,j<Length[allorder2higherCoef[[i]]]+1,j++,
		AppendTo[allCoefs,allorder2higherCoef[[i,j,2]]];
	];
];
Length[allCoefs]


(* ::Input::Initialization:: *)
(* Saving a_(s-1)i before evaluating their phi functions at 0 *)
asOriginal = a[[s-1]];
a20=a20/.L->0;
a21=a21/.L->0;
a22=a22/.L->0;
(* To append coefs from order3higher *)
(*   Applying conditions for order 4 holding weakly on last     *)
(*   intermediate stage:                                        *) 
(*                 a_(s-1)i(0) instead of a_(s-1)i(-hL)                 *)
(*  Always reread Theorem 4.7 in the Ho&Os05 paper when applying this step, "holding weakly" has a very particular meaning *) 
order3ahigher = {Sum[a[[s-1]][[i]]*c[[i]]^2/(2!),{i,stages-1}]==phi3[-L*h]};
order3bhigher={ExpandNonCommutativeMultiply[Sum[a[[s-1]][[i]]**J**(phi2c[[i]]*c[[i]]^2-Sum[a[[i-1]][[j]]*c[[j]],{j,2,i-1}]),{i,2,stages-1}]] == 0};
(* The following line takes care of psi_p(0)=0 in "holding weakly" *)
order3ahigher = order3ahigher/.{L->0};
allorder3higher = Join[order3ahigher,order3bhigher];
allorder3higherCoef = CoefficientRules[allorder3higher/.Equal->Subtract,termTypes];
For[i = 1, i<Length[allorder3higherCoef]+1,i++,
	For[j=1,j<Length[allorder3higherCoef[[i]]]+1,j++,
		AppendTo[allCoefs,allorder3higherCoef[[i,j,2]]];
	];
];
Length[allCoefs]


(* ::Input::Initialization:: *)
(* Don't run the following line unnesecarily because it takes a long time *)
(*solution=Solve[allCoefs==0]; *)
(* Importing the solution for the weights*)
solution = ToExpression[Import["/home/warhammer/Documents/Thoma_Zoto/UALBERTA_courses/Research/math/General_3_2_EXP_Method_half_three_fourth_one.txt", "List"
]];


(* ::Input::Initialization:: *)
(* Saving the 'solution' output for later use *)
(* Only run this if you are solving the system and not importing the solution.*)
(*Export["/home/warhammer/Documents/Thoma_Zoto/UALBERTA_courses/Research/math/General_3_2_EXP_Method_half_three_fourth_one.txt",solution];*)


(* ::Input::Initialization:: *)
(* Selecting the method *)
method=solution[[1]];
(*Variables[Values[method]];*)
a[[s-1]]=asOriginal;
(* Displaying the final Butcher tableau *)
a = a/.method;


(* ::Input::Initialization:: *)
(* Using the freedom in selecting coefficients *)
(* Making sure the second order method is never 3rd order *)
order3alower = {Sum[a[[s]][[i]]*c[[i]]^2/(2!),{i,stages}]!=phi3[-L*h]};
order3blower={Sum[a[[s]][[i]]*(phi2c[[i]]*c[[i]]^2-Sum[a[[i-1]][[j]]*c[[j]],{j,2,i-1}]),{i,2,stages}] != 0};
(* Stiff Conditions reduce to Classical Conditions when L=0, so we make sure that the 3nd order stiff method does NOT satisfy Classical Conditions of order 3. *)
order3alower=order3alower/.{L->0};
order3blower=order3blower/.{L->0};

allorder3lower = Join[order3alower, order3blower]/.method;

(* Alerts if one of the 4th order conditions is satisfied by the third order method without choosing the free variables*)
For[i=1,i<Length[allorder3lower]+1,i++,
	If[allorder3lower[[i]],,Print["Condition ", i, " is satisfied by second order methid. Check why."];];
];


(* ::Input::Initialization:: *)
(* This cell helps to minimize the the 4th order error coming from having to sattisfy the stiff conditions for order 3 strongly *)
order3ahigherO = {Sum[a[[s-1]][[i]]*c[[i]]^2/(2!),{i,stages-1}]==phi3[-L*h]};
order3bhigherO={ExpandNonCommutativeMultiply[Sum[a[[s-1]][[i]]**J**(phi2c[[i]]*c[[i]]^2-Sum[a[[i-1]][[j]]*c[[j]],{j,2,i-1}]),{i,2,stages-1}]] == 0};

allorder3higherO = Join[order3ahigherO,order3bhigherO];
(* Deletes conditions already sattisfied strongly*)
For[i=1,i<Length[allorder3higherO]+1,i++,
	If[allorder3higherO[[i]],allorder3higherO[[i]]=0;,Print["Condition ", i, " is False. Check why."];allorder3higherO=Delete[allorder3higherO,i];];
];

allorder3higherCoefO = CoefficientRules[allorder3higherO/.Equal->Subtract,termTypes];


(* ::Input::Initialization:: *)
(* This cell helps to minimize the the 4th order error coming from having to sattisfy the stiff conditions for order 4 weakly *)
(* Saving a_(s-2)i before evaluating their phi functions at 0 *)
(* a_(s-2)i will be used again for displaying the method at the end *)
asOriginal2 = a[[s-1]];
a[[s-1]][[1]]=a[[s-1]][[1]]/.L->0;
a[[s-1]][[2]]=a[[s-1]][[2]]/.L->0;
a[[s-1]][[3]]=a[[s-1]][[3]]/.L->0;
(*   Applying conditions for order 4 holding weakly on last     *)
(*   intermediate stage:                                        *) 
(*                 a_(s-1)i(0) instead of a_(s-1)i(-hL)                 *)
(*  Always reread Theorem 4.7 in the Ho&Os05 paper when applying this step, "holding weakly" has a very particular meaning *) 
order4ahigher = {Sum[a[[s-1]][[i]]*c[[i]]^3/(3!),{i,stages-1}]==phi4[-L*h]};
order4bhigher={ExpandNonCommutativeMultiply[Sum[a[[s-1]][[i]]**J**(phi3c[[i]]*c[[i]]^3-Sum[a[[i-1]][[j]]*c[[j]]^2/(2!),{j,2,i-1}]),{i,2,stages-1}]] == 0};
order4chigher = {ExpandNonCommutativeMultiply[Sum[a[[s-1]][[i]]**J**Sum[a[[i-1]][[j]]**J**(phi2c[[j]]*c[[j]]^2-Sum[a[[j-1]][[k]]*c[[k]],{k,2,j-1}]),{j,2,i-1}],{i,2,stages-1}]] == 0};
order4dhigher={ExpandNonCommutativeMultiply[Sum[a[[s-1]][[i]]*c[[i]]**K**(phi2c[[i]]*c[[i]]^2-Sum[a[[i-1]][[j]]*c[[j]],{j,2,i-1}]),{i,2,stages-1}]] == 0};
(* The following line takes care of psi_p(0)=0 in "holding weakly" *)
order4ahigher = order4ahigher/.{L->0};
allorder4higher = Join[order4ahigher,order4bhigher,order4chigher,order4dhigher];

(* Deletes conditions already sattisfied strongly *)
For[i=1,i<Length[allorder4higher]+1,i++,
	If[allorder4higher[[i]],allorder4higher[[i]]=0;,Print["Condition ", i, " is False. Check why."];allorder4higher=Delete[allorder4higher,i];];
];

allorder4higherCoef = CoefficientRules[allorder4higher/.Equal->Subtract,termTypes];


(* ::Input::Initialization:: *)
(* We will put the coefficients in front of all term types involving J,K and W. normSquared is the squared 2-norm of that vector *) 
normSquared = 0;
For[i = 1, i<Length[allorder3higherCoefO]+1,i++,
	For[j=1,j<Length[allorder3higherCoefO[[i]]]+1,j++,
	normSquared = normSquared + (allorder3higherCoefO[[i,j,2]])^2;
	]; 
	];

For[i = 1, i<Length[allorder4higherCoef]+1,i++,
	For[j=1,j<Length[allorder4higherCoef[[i]]]+1,j++,

normSquared = normSquared + (allorder4higherCoef[[i,j,2]])^2;
	];
];



(* ::Input::Initialization:: *)
(*Minimizing the norm of the coefficients in front of terms involving J,K,W in the stiff order conditions for order 5 *)
norm = {Sqrt[normSquared]};
optimalParameters = NMinimize[Join[norm,allorder3lower],Variables[Values[method]]]
(* Approximating the selected parameters with fractions *)
fractionParams=Rationalize[optimalParameters[[2]],0.1]


(* ::Input::Initialization:: *)
(* Displaying the final Butcher tableau *)
(* normSquared for the method with fractional coefficients is approximately 4.81 *)
a[[s-1]] = asOriginal2;
aERK32ZBs = a/.fractionParams
aERK32ZBs[[4]][[1]];
aERK32ZBs[[4]][[2]];
aERK32ZBs[[4]][[3]];
aERK32ZBs[[4]][[4]];


(* ::Input:: *)
(*aERK32ZBs/.{L->0}*)
(**)


(* ::Input::Initialization:: *)
testCell = NotebookRead[Cells[][[39]]];
CellToTeX[testCell]

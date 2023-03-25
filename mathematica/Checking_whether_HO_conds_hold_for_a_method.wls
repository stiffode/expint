#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input::Initialization:: *)
(*NCAlgebra allows for non commutative expr*)
PacletInstall["https://github.com/NCAlgebra/NC/blob/v6.0.0/NCAlgebra-6.0.0.paclet?raw=true"];
<<NCAlgebra`;



(* ::Input::Initialization:: *)
(* To clear all variables *)
Clear["Global`*"]

(*Modify the number of stages*)
stages = 4
s=stages;
(*Building the Butcher Tableau and vector of c_i*)
a=Table[Symbol["a"<>ToString[i]<>ToString[j]],{i,0,stages-1},{j,0,i}]
c=Table[Symbol["c"<>ToString[i]],{i,0,stages-1}]
(*Evaluating each phi function at each c_i*)
phi1c=Table[phi1[-c[[i]]*L*h],{i,1,stages}];
phi2c=Table[phi2[-c[[i]]*L*h],{i,1,stages}];
phi3c=Table[phi3[-c[[i]]*L*h],{i,1,stages}];
phi4c=Table[phi4[-c[[i]]*L*h],{i,1,stages}];

(*Stiff conditions as in 05 Hochbruck paper*)SetNonCommutative[J,K,L,W,a]
order1={Sum[a[[s]][[i]],{i,stages}]==phi1[-L*h]};
order2a={Sum[a[[s]][[i]]*c[[i]],{i,stages}]==phi2[-L*h]};
order2b=Table[Sum[a[[i]][[k]],{k,1,i}]==phi1c[[i+1]]*c[[i+1]],{i,1,stages-1}];
order3a={Sum[a[[s]][[i]]*c[[i]]^2/(2!),{i,1,stages}]==phi3[-L*h]};
order3b={Sum[a[[s]][[i]]**J**(phi2c[[i]]*c[[i]]^2-Sum[a[[i-1]][[j]]*c[[j]],{j,2,i-1}]),{i,2,stages}]==0};
order4a={Sum[a[[s]][[i]]*c[[i]]^3/(3!),{i,stages}]==phi4[-L*h]};
order4b={Sum[a[[s]][[i]]**J**(phi3c[[i]]*c[[i]]^3-Sum[a[[i-1]][[j]]*c[[j]]^2/(2!),{j,2,i-1}]),{i,2,stages}]==0};
order4c={Sum[a[[s]][[i]]**J**Sum[a[[i-1]][[j]]**J**(phi2c[[j]]*c[[j]]^2-Sum[a[[j-1]][[k]]*c[[k]],{k,2,j-1}]),{j,2,i-1}],{i,2,stages}]==0};
order4d={Sum[a[[s]][[i]]*c[[i]]**K**(phi2c[[i]]*c[[i]]^2-Sum[a[[i-1]][[j]]*c[[j]],{j,2,i-1}]),{i,2,stages}]==0};

(* Stiff order conditions from Luan14 paper*)
order5a = {Sum[a[[s]][[i]]*c[[i]]^4/(4!),{i,stages}]==phi5[-L*h]};
order5b={Sum[a[[s]][[i]]**J**(phi4c[[i]]*c[[i]]^4-Sum[a[[i-1]][[j]]*c[[j]]^3/(3!),{j,2,i-1}]),{i,2,stages}] == 0};
order5c = {Sum[a[[s]][[i]]**J**Sum[a[[i-1]][[j]]**J**(phi3c[[j]]*c[[j]]^3-Sum[a[[j-1]][[k]]*c[[k]]^2/(2!),{k,2,j-1}]),{j,2,i-1}],{i,2,stages}] == 0};
order5d = {Sum[a[[s]][[i]]**J**Sum[a[[i-1]][[j]]**J**Sum[a[[j-1]][[k]]**J**(phi2c[[k]]*c[[k]]^2-Sum[a[[k-1]][[l]]*c[[l]],{l,2,k-1}]),{k,2,j-1}],{j,2,i-1}],{i,2,stages}] == 0};
order5e = {Sum[a[[s]][[i]]**J**Sum[(a[[i-1]][[j]]*c[[j]])**K**((phi2c[[j]]*c[[j]]^2 )-Sum[a[[j-1]][[k]]*c[[k]],{k,2,j-1}]),{j,2,i-1}],{i,2,stages}]== 0};
order5f={Sum[a[[s]][[i]]*c[[i]]**K**(phi3c[[i]]*c[[i]]^3-Sum[a[[i-1]][[j]]*c[[j]]^2/(2!),{j,2,i-1}]),{i,2,stages}] == 0};
order5g = {Sum[a[[s]][[i]]*c[[i]]**K**Sum[a[[i-1]][[j]]**J**(phi2c[[j]]*c[[j]]^2-Sum[a[[j-1]][[k]]*c[[k]],{k,2,j-1}]),{j,2,i-1}],{i,2,stages}] == 0};

(* Condition 17 is satisfied like it is satisfied in the paper by Lu&Os14 *)
order5h=Table[a[[s]][[i]]==0 ||phi2c[[i]]*c[[i]]^2-Sum[a[[i-1]][[j]]*c[[j]],{j,2,i-1}]==0,{i,2,stages}] ;

order5i={Sum[a[[s]][[i]]*c[[i]]^2**W**(phi2c[[i]]*c[[i]]^2-Sum[a[[i-1]][[j]]*c[[j]],{j,2,i-1}]),{i,2,stages}] == 0};


(* ::Input::Initialization:: *)
(* The third order method of ERK43ZB *)
c0=0;c1=1/6;c2=1/2;c3=1/2;c4=1;
a00=1/6*phi1[-L*h/6];

a11=3/2*phi2[-L*h/2]+1/2*phi2[-L*h/6];
a10=1/2*phi1[-L*h/2] -a11;

a22 = -19/180*phi1[-L*h]-1/6*phi1[-L*h/2]-1/6*phi1[-L*h/6]-1/6*phi2[-L*h/2]+1/9*phi2[-L*h/6]-1/5*phi3[-L*h/2];
a21 =  19/60*phi1[-L*h]+1/2*phi1[-L*h/2]+1/2*phi1[-L*h/6]+2*phi2[-L*h/2]+13/6*phi2[-L*h/6]+3/5*phi3[-L*h/2];
a20 = 1/2*phi1[-L*h/2] - a21-a22;

a33=phi2[-L*h]+phi2[-L*h/2]-6*phi3[-L*h]-3*phi3[-L*h/2];
a31=3*phi2[-L*h]-9/2*phi2[-L*h/2]-5/2*phi2[-L*h/6]+6*a33+a21;
a32=6*phi3[-L*h]+3*phi3[-L*h/2]-2*a33+a22;
a30=phi1[-L*h]-a31-a32-a33;

(*  Checking whether conditions hold exactly for orders 1,2  *)
allorder1 = order1;
FullSimplify[allorder1]
allorder2 = Join[order2a,order2b];
FullSimplify[allorder2]

(*Checking whether condition 3 holds weakly i.e.b_i(0) instead of b_i(-hL)*)(*Always reread Theorem 4.7 in the 05 paper when applying this step,"holding weakly" has a very particular meaning*)a30=a30/. L->0;
a31=a31/. L->0;
a32=a32/. L->0;
a33=a33/. L->0;
phi1[0]=1;phi2[0]=1/2;phi3[0]=1/(3!);phi4[0]=1/(4!);
phi5[0]=1/(5!);

(*The following line takes care of psi_p(0)=0 in "holding weakly"*)

order3a=order3a/. {L->0};
allorder3=Join[order3a,order3b];

FullSimplify[allorder3]

(* Checking that this method is never 4th order *)
allorder4=Join[order4a,order4b,order4c,order4d];
allorder4/.{L->0}




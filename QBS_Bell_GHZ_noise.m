close all;clear;clc;

L0=[1;0];
L1=[0;1];

L00=kron(L0,L0);
L01=kron(L0,L1);
L10=kron(L1,L0);
L11=kron(L1,L1);

L000=kron(L0,L00);
L001=kron(L0,L01);
L010=kron(L0,L10);
L011=kron(L0,L11);
L100=kron(L1,L00);
L101=kron(L1,L01);
L110=kron(L1,L10);
L111=kron(L1,L11);

%initialization
GHZ=sqrt(2)*(L000+L111);
Bell00=sqrt(2)*(L00+L11);
rhoi=kron(L00*L00',kron())

%quantum gate


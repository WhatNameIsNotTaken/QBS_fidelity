function myKron=mykron(a,b,c,d,e,f,g)
if nargin==1
    myKron=kron(a,a);

elseif nargin==2
  myKron=kron(a,b);

elseif nargin==3
 myKron=kron(kron(a,b),c);

 elseif nargin==4
 myKron=kron(kron(kron(a,b),c),d);
 
 elseif nargin==5
 myKron=kron(kron(kron(kron(a,b),c),d),e);
 
 elseif nargin==6
 myKron=kron(kron(kron(kron(kron(a,b),c),d),e),f);
 
  elseif nargin==7
 myKron=kron(kron(kron(kron(kron(kron(a,b),c),d),e),f),g);
end
global L
global S
L=30 //Taille de la boite
S=L*L


//Laplacien en dimension 2


//parametrage pour pouvoir parcourrir les éléments de la matrice sur une surface
function y=code12(i,j)
    y=L*(i-1)+j
endfunction

function [i,j]=code21(y)
    i=modulo(y,L)
    j=y/L
endfunction


//construction de la matrice.
H=zeros(S,S)

for i=1:L
    for j=1:L
        i2=modulo(i,L)+1
        j2=modulo(j,L)+1
        
        i3=modulo(i-2+L,L)+1
        j3=modulo(j-2+L,L)+1
        H(code12(i,j),code12(i,j2))=1
        H(code12(i,j),code12(i,j3))=1
        H(code12(i,j),code12(i2,j))=1
        H(code12(i,j),code12(i3,j))=1
    end
end

//on diagonalise
[D,P]=bdiag(H)

//on affiche un vecteur propre

x=P(:,6)
xPlot=matrix(x,L,L)

clf()
plot3d(1:L,1:L,100*xPlot)
a=gca()
xs2png(gcf(),'onde_delocalise_sans_desordre.png');
//Sa.data_bounds=[0,0,min(x);L,L,max(x)]

global L
global S
global w
L=30 //Taille de la boite
S=L*L //la surface de la boite en dimension 2
w=5 //force du desordre.

//il faut apparament augmenter la mémoire disponible.
//stacksize(100000000)


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

//on ajoute l'aléa.
alea=rand(1,S)
W=w*diag(alea)

H=H+W

//on diagonalise
[D,P]=bdiag(H)

//On calcule le répartition des valeurs propres
dd=diag(D)
p=zeros(1,20)
eX=zeros(1,20)
for u=1:20
    
    a=find(dd>(u*(w+6)/(20-2)-4))
    d2=dd(a)
    b=find(d2<((u+1)*(w+6)/(20-2)-4))
    p(u)=length(b)
    eX(u)=floor(10*(u*(w+6)/(20-2)-4))/10
end

//on affiche un vecteur propre

x=P(:,6)
xPlot=matrix(x,L,L)
if(%f)
    clf()
    plot3d(1:L,1:L,100*xPlot)
    a=gca()
end

//On sauvegarde l'image.
//xs2png(gcf(),'onde_delocalise_sans_desordre.png');
//Sa.data_bounds=[0,0,min(x);L,L,max(x)]

//On coupe une sous boite:
Hsb=H
k=5
compt=0
for compteur1=1:k
    for compteur2=1:k
        limiteUpI=L/k*compteur1
        limiteDownI=L/k*(compteur1-1)+1
        limiteUpJ=L/k*(compteur2)
        limiteDownJ=L/k*(compteur2-1)+1
        
        for i=(limiteDownI):(limiteUpI)
        for j=(limiteDownJ):(limiteUpJ)
        i2=modulo(i,L)+1
        j2=modulo(j,L)+1
        
        i3=modulo(i-2+L,L)+1
        j3=modulo(j-2+L,L)+1
        if(i2>limiteUpI)
        Hsb(code12(i,j),code12(i2,j))=0
        Hsb(code12(i2,j),code12(i,j))=0
       compt=compt+1
        end
        if(i3<limiteDownI)
            Hsb(code12(i,j),code12(i3,j))=0
            Hsb(code12(i3,j),code12(i,j))=0
            compt=compt+1
        end
        if(j2>limiteUpJ)
        Hsb(code12(i,j),code12(i,j2))=0
        Hsb(code12(i,j2),code12(i,j))=0
        compt=compt+1
        end
        if(j3<limiteDownJ )
            Hsb(code12(i,j),code12(i,j3))=0
            Hsb(code12(i,j3),code12(i,j))=0
            compt=compt+1
        end
        end
    end
    end
    
    
end




//On choisit une valeur propre.
vp=1
//On désigne un vecteur test
pointA=zeros(S,1)
pointA(S/2+L/2,1)=1

pointAPlot=matrix(pointA,L,L)

//Puis on calcule la fonction de green à partir de ce point.
GreenVp=(Hsb-vp*eye(S,S))^(-1)
GreenA=GreenVp*pointA
//Et on l'affiche
greenAPlot=matrix(GreenA,L,L)
figure(2)
clf()
plot3d(1:L,1:L,10*greenAPlot)

//Puis on calcule la fonction de green à partir de ce point.
GreenVpH=(H-vp*eye(S,S))^(-1)
GreenAH=GreenVpH*pointA
//Et on l'affiche
greenAHPlot=matrix(GreenAH,L,L)
figure(3)
clf()
plot3d(1:L,1:L,10*greenAHPlot)


if(%f)  //Pour éviter de réaliser cette partie.

//Et on l'affiche le potentiel aléatoire
//WPlot=matrix(alea,L,L)
WPlot=zeros(L,L)
figure(5)
clf()


//on affiche en rouge
f=gcf()
xset("default")
plot3d(1:L,1:L,10*WPlot)

//on affiche la densité d'état

figure(6)
clf()
bar(eX,p)

end

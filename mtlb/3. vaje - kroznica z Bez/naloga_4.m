addpath('..\1. vaje - bernsteinova baza, deCasteljau,Bezier', '..\2. vaje');

% 4. naloga:

n = 5;
fi = pi/n;

% izračunamo kontrolni poligon na vse načine:
B1 = bezierarc(fi,1);
B2 = bezierarc(fi,2);
B3 = bezierarc(fi,3);

% izračun točk na Bezierjevih krivuljah
t = linspace(0,1,1001);
b1 = bezier(B1,t);
b2 = bezier(B2,t);
b3 = bezier(B3,t);

% radialne napake max(|1-b(t)|, t \in [0,1])
[e1,e2,e3] = deal(0); % vrednost vseh treh spremenljivk nastavi na 0
for i = 1:length(t)
    e1 = max(e1,abs(1-norm(b1(i,:))));
    e2 = max(e2,abs(1-norm(b2(i,:))));
    e3 = max(e3,abs(1-norm(b3(i,:))));
end
disp('napake = ')
[e1,e2,e3]

% rotacijska matrika
c = cos(2*fi);
s = sin(2*fi);
R = [c -s; s c];

% narišemo krožnice:
clf
hold on
plot(cos(2*pi*t),sin(2*pi*t),'Color',0.7*[1 1 1],'LineWidth',4);
plot(b1(:,1),b1(:,2),'r');
plot(b2(:,1),b2(:,2),'g');
plot(b3(:,1),b3(:,2),'b');
for i = 1:n-1
    B1 = (R*B1')';
    B2 = (R*B2')';
    B3 = (R*B3')';
    b1 = bezier(B1,t);
    b2 = bezier(B2,t);
    b3 = bezier(B3,t);
    plot(b1(:,1),b1(:,2),'r');
    plot(b2(:,1),b2(:,2),'g');
    plot(b3(:,1),b3(:,2),'b');
end
hold off
axis equal

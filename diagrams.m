 clc
 clear all
 close all

fileID = fopen('data.txt', 'r');
sizes = zeros(1,4);

% B
sizes(1) = fscanf(fileID, '%d\n', 1);
Jacobi_B = fscanf(fileID, '%f\n', sizes(1));
sizes(2) = fscanf(fileID, '%d\n', 1);
GS_B = fscanf(fileID, '%f\n', sizes(2));

% C
sizes(3) = fscanf(fileID, '%d\n', 1);
Jacobi_C = fscanf(fileID, '%f\n', sizes(3));
sizes(4) = fscanf(fileID, '%d\n', 1);
GS_C = fscanf(fileID, '%f\n', sizes(4));

f1 = figure;
semilogy(1:length(Jacobi_B), Jacobi_B)
hold on
semilogy(1:length(GS_B), GS_B)
legend('metoda Jacobiego', 'metoda Gaussa-Seidla')
xlabel('Nr iteracji')
ylabel('Wartość normy z residuum')
print("normaB", "-dpng")
hold off

f2 = figure;
semilogy(1:length(Jacobi_C), Jacobi_C)
hold on
semilogy(1:length(GS_C), GS_C)
legend('metoda Jacobiego', 'metoda Gaussa-Seidla')
xlabel('Nr iteracji')
ylabel('Wartość normy z residuum')
print("normaC", "-dpng")
hold off


% E
i = 1;
sizes = zeros(1, 6);
durationsJacobi = zeros(1, 6);
durationsGS = zeros(1, 6);
durationsLU = zeros(1, 6);

while(i < 7)
    sizes(i) = fscanf(fileID, '%d\n', 1);
    durationsJacobi(i) = fscanf(fileID, '%f\n', 1);
    durationsGS(i) = fscanf(fileID, '%f\n', 1);
    durationsLU(i) = fscanf(fileID, '%f\n', 1);
    i = i + 1;
end
 
fclose(fileID);

f3 = figure;
plot(sizes, durationsJacobi, '-o')
hold on
grid on
plot(sizes, durationsGS, '-o')
plot(sizes, durationsLU, '-o')
xlabel('Liczba niewiadomych')
ylabel('Czas trwania [s]')
legend('metoda Jacobiego', 'metoda Gaussa-Seidla', 'metoda faktoryzacji LU', 'Location', 'northwest')
print("zadE", "-dpng")

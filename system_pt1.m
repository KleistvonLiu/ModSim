% MODSIM Laborpraktikum, 1. Aufgabe
%
% Dr.-Ing. Th. Range, Dr.-Ing. S. Dyblenko
%
% zu erg?nzende Codezeilen sind mit ">>> erg?nzen ...." gekennzeichnet
% Berechnung des Systems "PT1-Glied"
%
% (Hinweis: Die Struktur des Programms erlaubt eine Einbindung in
% SIMULINK als sog. S-Function; dies ist jedoch im Rahmen des
% Praktikums nicht vorgesehen!)
function [sys, x0] = system_pt1 ( t, x, u, flag )
% Eingabe-Parameter:
% t - Zeit
% x - x(t) Zustandsvektor zum Zeitpunkt t
% u - u(t) Vektor der Eingangssignale zum Zeitpnukt t
% flag - Steuerparameter, legt fest, welche Ausgabe
% gefordert wird:
Tm = 10; % Zeitkonstante des PT1-Gliedes
if flag == 0 % Ausgabe der Anfangswerte f¨¹r den
        % Zustand auf Vektor x0
    x0 = 0; 
    sys = [1,0,1,1,0,0];% diese Zeile ist nur f¨¹r Simulink n?tig,
        % sie gilt so NUR in diesem Beispiel!
elseif abs(flag) == 1 % Ausgabe der Ableitungen von x auf
        % Vektor sys = x' = f(x(t),u(t),t)
    sys = -1/Tm*x+1/Tm*u;
elseif flag == 3 % Ausgabe der Ausgangswerte des
        % Systems auf Vektor sys:
        % sys = y = g(x(t),u(t),t)
    sys = u*(1-exp(-t/Tm));
else % bei anderen Flagwerten nichts ausgeben
    sys = [];
end
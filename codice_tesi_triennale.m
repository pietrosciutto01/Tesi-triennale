%% CODICE MATLAB PER LA TESI TRIENNALE - PIETRO SCIUTTO (5017366)
% SVILUPPO DI UN CODICE DI CALCOLO PER IL DIMENSIONAMENTO DI UNA RETE DI
% DRENAGGIO URBANO CON IL METODO DEL VOLUME DI INVASO

clear; clc; close all

format bank % per avere 2 cifre decimali dopo la virgola nella command window
% format long % per avere 15 cifre decimali dopo la virgola nella command window

numero_figure = 0;


%% dati di partenza

w0 = 40; % contributi superficiali - volume dei piccoli invasi (m3 / ha)
a = 70; % parametro 'a(T)' della LSPP (mm/h^n)
n = 0.30; % parametro 'n' della LSPP (-)
ks = 70; % coefficiente di scabrezza di gauckler-strickler delle condotte (m1/3 s-1)
diametri = [0.3 0.4 0.5 0.6 0.7 0.8 0.9]; % diametri delle condotte


%% disp dei dati di partenza

disp('-----------------------------------------------------------------')
disp('DATI DI PARTENZA:')
disp([' • il contributo superficiale - piccoli invasi è w0 = ' num2str(w0) ' m3/ha'])
disp([' • il parametro ''a'' della LSPP è a = ' num2str(a) ' mm/h^n'])
disp([' • il parametro ''n'' della LSPP è a = ' num2str(n) ' (-)'])
disp([' • il coefficiente di scabrezza di gauckler-strickler delle condotte è ks = ' num2str(ks) ' m^(1/3)/s^(-1)'])
disp([' • i diametri commerciali delle condotte disponibili per il progetto sono: ' sprintf('%.2f ', diametri) 'metri'])


%% importazione dei dati da QGIS

% lettura dello shapefile con la planimetria definitiva della rete
S_completa = shaperead("rete_ac_def_dissolto_generalized.shp");
campi_da_estrarre = {'Fields', 'X', 'Y', 'fid'}; % estraggo dallo shapefile solo le colonne che mi interessano (fid, X, Y)
S = rmfield(S_completa, setdiff(fieldnames(S_completa), campi_da_estrarre)); % creo una nuova struct uguale a quella di partenza ma solo con i dati utili
% emilinazione dei 'NaN' in fondo ai vettori delle coordinate X e Y
for i = 1 : length(S)
    S(i).X = S(i).X(1:end-1);
    S(i).Y = S(i).Y(1:end-1);
end

% importazione del DTM
[DTM, info_DTM] = readgeoraster('DTM_merged_filled copia.tif');
x0_DTM = info_DTM.XWorldLimits(1); y0_DTM = info_DTM.YWorldLimits(2);
dX_DTM = info_DTM.CellExtentInWorldX; dY_DTM = -(info_DTM.CellExtentInWorldY); % lati della cella quadrata di lato 1m

% importazione della mappa di accumulazione
[accumulazione, info_accumulazione] = readgeoraster('accumulazione_r_flow_ritagliata.tif'); % per determinare le aree a monte dei rami sorgente

% importazione della carta di uso del suolo
[uso_suolo, info_uso_suolo] = readgeoraster('uso_suolo_raster_ritagliato.tif');
x0_uso_suolo = info_uso_suolo.XWorldLimits(1); y0_uso_suolo = info_uso_suolo.YWorldLimits(2);
dX_uso_suolo = info_uso_suolo.CellExtentInWorldX; dY_uso_suolo = -(info_uso_suolo.CellExtentInWorldY);

% importazione dell'ortofoto
[orto_merged, R_merged] = readgeoraster('prov2.tif');
% se è in scala di grigi, l'ortofoto viene replicata su 3 canali per farla diventare RGB
if size(orto_merged, 3) == 1
    orto_merged = repmat(orto_merged, [1 1 3]); % conversione da scala di grigi a RGB
elseif size(orto_merged, 3) > 3
    orto_merged = orto_merged(:, :, 1:3); % si mantengono solo i primi 3 canali
end


%% vettore contenente i rami terminali

rami_terminali = [19 33 55 104 207]; % questo vettore deve essere scritto MANUALMENTE


%% associo a tutti i punti (X,Y) la quota z

% x = x0 + i * dX → i = (x - x0)/dX | y = y0 + i * dY → i = (y - y0)/dY
% 'round' → arrotondare al decimale o al numero intero più vicino

for i = 1:length(S)
    for j = 1:length(S(i).X)
        ix = round((S(i).X(j) - x0_DTM) / dX_DTM);
        iy = round((S(i).Y(j) - y0_DTM) / dY_DTM);
        S(i).Z(j) = DTM(iy, ix);
    end
end


%% calcolo della lunghezza 3D dei rami

for i = 1:length(S)
    S(i).L = 0;
    for j = 1:length(S(i).X)
        if j < length(S(i).X)
            % lunghezza del ramo: L = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2) (m)
            S(i).L = S(i).L + sqrt(((S(i).X(j+1) - S(i).X(j))^2 + (S(i).Y(j+1) - S(i).Y(j))^2 + (S(i).Z(j+1) - S(i).Z(j))^2));
        end
    end
end


%% eliminazione delle righe corrispondenti ai rami con L < 1m
i = 1;
while i < length(S)
    if S(i).L < 1
        S(i,:) = [];
    else
        i = i + 1;
    end
end


%% disp delle statistiche sulle lunghezze L dei rami

disp('-----------------------------------------------------------------')
disp('LUNGHEZZE L DEI RAMI:')
valori_L = [S.L];
valori_fid = [S.fid];
[~, indice_max] = max(valori_L); % indice del valore massimo nella colonna L
max_fid = valori_fid(indice_max); % valore corrispondente nella colonna fid
max_L = valori_L(indice_max); % valore massimo nella colonna L
disp([' • il ramo ' num2str(max_fid) ' è il ramo più lungo della rete, con una lunghezza L = ' num2str(max_L,'%.2f') ' m']);

L_media = mean(valori_L);
disp([' • la lunghezza media dei rami della rete è pari a ' num2str(L_media,'%.2f') ' m'])

[~, indice_min] = min(valori_L);
min_fid = valori_fid(indice_min);
min_L = valori_L(indice_min);
disp([' • ll ramo ' num2str(min_fid) ' è il ramo più corto della rete, con una lunghezza L = ' num2str(min_L,'%.2f') ' m']);


%% determinazione dei punti di inizio e fine di ogni ramo

% si parte dai rami terminali e si risale la rete, fino ad arrivare ai rami sorgente
rami_visitati = []; % serve per tenere traccia dei rami già visitati per evitare di ricontrollarli
tolleranza_approx = 10^(-3); % il punto di fine ramo e quello di inizio ramo del successivo possono non coincidere perfettamente, ma...
% essere estremamente vicini: per questo si assume una certa tolleranza in questo passaggio

for i = 1:length(rami_terminali)

    % vengono assegnate le coordinate del punto iniziale e finale dei rami terminali
    % facendo l'ipotesi che il primo punto sia l'inizio e il punto
    % end la fine; successivamente l'ipotesi verrà verificata, e se dovesse
    % riverlarsi errata verrebbe invertito l'ordine del ramo inizialmente ipotizzato

    indice_ramo = find([S.fid] == rami_terminali(i));
    S(indice_ramo).Inizio = [S(indice_ramo).X(1), S(indice_ramo).Y(1)];
    S(indice_ramo).Fine = [S(indice_ramo).X(end), S(indice_ramo).Y(end)];

    S(indice_ramo).convergein = 0;

    S(indice_ramo).id_sottobacino = char('A' + (i - 1)); % per la suddivisione in sottobacini

    for j = 1:length(S)
        if ((S(indice_ramo).X(end) == S(j).X(1) && S(indice_ramo).Y(end) == S(j).Y(1)) ||...
                (S(indice_ramo).X(end) == S(j).X(end) && S(indice_ramo).Y(end) == S(j).Y(end))) && ...
                S(j).fid ~= rami_terminali(i)

            % se il punto che è stato ipotizzato come fine invece è in
            % comune con almeno un altro ramo, allora viene invertito l'ordine
            % ipotizzato precedentemente
            S(indice_ramo).Inizio = [S(indice_ramo).X(end), S(indice_ramo).Y(end)];
            S(indice_ramo).Fine = [S(indice_ramo).X(1), S(indice_ramo).Y(1)];
        end
    end

    % determinazione di inizio e fine dei rami che non siano terminali:

    rami_adiacenti_i = []; % vettore di appoggio per inserire i rami adiacenti a quello che si sta controllando
    rami_adiacenti_i = [rami_adiacenti_i; rami_terminali(i)];

    rami_visitati = [rami_visitati; rami_terminali(i)];

    while ~isempty(rami_adiacenti_i)
        tmp = 0;
        for j = 1:length(S)
            if sum(rami_visitati == S(j).fid) == 0 % si evita di controllare il ramo attuale con se stesso & non lo si considera se è già stato visitato prima

                indice_ramo = find([S.fid] == rami_adiacenti_i(1));

                if (abs((S(j).X(1) - S(indice_ramo).Inizio(1))) < tolleranza_approx) && (abs(S(j).Y(1) - S(indice_ramo).Inizio(2)) < tolleranza_approx)
                    rami_adiacenti_i = [rami_adiacenti_i; S(j).fid];
                    S(j).Fine = [S(j).X(1),S(j).Y(1)];
                    S(j).Inizio = [S(j).X(end),S(j).Y(end)];

                    S(j).X = [S(j).X(end:-1:1)];
                    S(j).Y = [S(j).Y(end:-1:1)];
                    S(j).Z = [S(j).Z(end:-1:1)];

                    S(j).convergein = rami_adiacenti_i(1); % si scrive in questa colonna il ramo in cui converge il ramo j-esimo

                    S(j).id_sottobacino = S(indice_ramo).id_sottobacino;

                    rami_visitati = [rami_visitati; S(j).fid];
                    tmp = tmp + 1;
                elseif (abs(S(j).X(end) - S(indice_ramo).Inizio(1)) < tolleranza_approx) && (abs(S(j).Y(end) - S(indice_ramo).Inizio(2)) < tolleranza_approx)
                    rami_adiacenti_i = [rami_adiacenti_i; S(j).fid];
                    S(j).Fine = [S(j).X(end),S(j).Y(end)];
                    S(j).Inizio = [S(j).X(1),S(j).Y(1)];

                    S(j).convergein = rami_adiacenti_i(1);

                    S(j).id_sottobacino = S(indice_ramo).id_sottobacino;

                    rami_visitati = [rami_visitati; S(j).fid];
                    tmp = tmp + 1;
                end
            end
        end

        indice_ramo = find([S.fid] == rami_adiacenti_i(1));

        if tmp == 0
            if sum(rami_terminali == rami_adiacenti_i(1)) == 0
                S(indice_ramo).sorgenteOcollettore = "sorgente";
            end
        else
            S(indice_ramo).sorgenteOcollettore = "collettore";
        end
        rami_adiacenti_i(1) = []; % viene rimosso il ramo appena visitato
    end
    indice_ramo = find([S.fid] == rami_terminali(i));
    S(indice_ramo).sorgenteOcollettore = 'terminale';

end


%% controllo per vedere che tutti i rami siano classificati

contatore = 0;
for i = 1 : length(S)
    if isempty(S(i).sorgenteOcollettore)
        contatore = contatore + 1;
    end
end

disp('-----------------------------------------------------------------')
disp('CLASSIFICAZIONE RAMI:')
if contatore == 0
    disp('tutti i rami sono classificati come S o C')
else
    disp([num2str(contatore) ' rami NON sono classificati come S o C'])
end


%% determinazione del n° di strahler di ogni ramo

rami_da_visitare = [S.fid]';
i = 1;
n_rami_rimanenti = length(S);
while ~isempty(rami_da_visitare)
    indice_ramo = find([S.fid] == rami_da_visitare(i));
    if S(indice_ramo).sorgenteOcollettore == "sorgente"
        S(indice_ramo).n_strahler = 1; % i rami sorgente hanno n_strahler = 1
        rami_da_visitare(i) = [];
        n_rami_rimanenti = n_rami_rimanenti - 1;
    else

        convergenti_in_i = [find([S.convergein] == S(indice_ramo).fid)];

        stop = 0; % variabile che serve per capire quando si sta considerando un ramo di cui non conosco l'area a monte

        strahler_tmp = []; % creazione di un vettore vuoto in cui mettere tutti i numeri di strahler dei rami convergenti in k
        for j = 1:length(convergenti_in_i)

            if isempty(S(convergenti_in_i(j)).n_strahler) % se il ramo non ha ancora n_strahler, allora stop = 1 e si esce dal for
                stop = 1;
                break
            else
                strahler_tmp = [strahler_tmp, S(convergenti_in_i(j)).n_strahler]; % altrimenti si aggiunge n_strahler del ramo considerato al vettore strahler_tmp
            end
        end
        % una volta uscito dal for, se tutti i rami convergenti in quello
        % considerato avevano l'area a monte e numero di strahler già calcolati, allora stop = 0
        % e si può calcolare l'area monte e il numero di strahler del ramo considerato; altrimenti
        % si passa al ramo successivo incrementando k

        if stop == 0
            strahler_tutti_uguali = 0;
            for j = 2 : length(strahler_tmp)
                if strahler_tmp(j) ~= strahler_tmp(1)
                    strahler_tutti_uguali = 1;
                end
            end
            % determinazione del numero di strahler:
            if strahler_tutti_uguali == 0
                S(indice_ramo).n_strahler = strahler_tmp(1) + 1;
            else
                S(indice_ramo).n_strahler = max(strahler_tmp);
            end

            rami_da_visitare(i) = [];
            n_rami_rimanenti = n_rami_rimanenti - 1;
        else
            i = i + 1;
        end
    end
    % se si è arrivati alla fine dei rami, allora si mette k = 1 e si
    % ricomincia da capo per calcolare le aree dei rami non calcolati precedentemente;
    % se non ci sono più rami da visitare, allora esce dal while
    if i == n_rami_rimanenti + 1
        i = 1;
    end
end


%% coefficiente di afflusso

coeff_afflusso = table2array(readtable('coeff_afflusso.xlsx'));

for i = 1:length(S)
    x_min(i) = round(min(S(i).X) - x0_uso_suolo) / dX_uso_suolo;
    x_max(i) = round(max(S(i).X) - x0_uso_suolo) / dX_uso_suolo;
    y_min(i) = round(min(S(i).Y) - y0_uso_suolo) / dY_uso_suolo;
    y_max(i) = round(max(S(i).Y) - y0_uso_suolo) / dY_uso_suolo;
    sottomatrice = uso_suolo(y_max(i) : y_min(i), x_min(i) : x_max(i));

    somma_valori = 0;

    for j = 1:size(sottomatrice, 1)
        for k = 1:size(sottomatrice, 2)

            num_cifre = floor(log10(sottomatrice(j,k))) + 1; % si calcola il numero di cifre nel numero
            primiDueCaratteri = floor(sottomatrice(j,k) / 10^(num_cifre - 2)); % si estraggono le due prime cifre

            valore_coeff_afflusso = 0;
            for h = 1 : size(coeff_afflusso, 1)
                if coeff_afflusso(h,1) == primiDueCaratteri
                    valore_coeff_afflusso = coeff_afflusso(h,2);
                    break
                end
            end

            somma_valori = somma_valori + valore_coeff_afflusso;
        end
    end
    S(i).fi = somma_valori / (size(sottomatrice,1) * size(sottomatrice,2));
end


%% calcolo dell'area dei rami

% area totale del bacino di interesse (ricavata da QGIS):
A_tot_mappa = 2812539 / (10^4); % (da m2 a ha)

% calcolo dell'area a monte dei rami sorgente
for i = 1:length(S)

    % si associa la quota z (m) ai punti di inizio e fine
    i_x_inizio = round((S(i).Inizio(1) - x0_DTM)/dX_DTM);
    i_y_inizio = round((S(i).Inizio(2) - y0_DTM)/dY_DTM);
    S(i).Inizio(3) = DTM(i_y_inizio,i_x_inizio);

    i_x_fine = round((S(i).Fine(1) - x0_DTM)/dX_DTM);
    i_y_fine = round((S(i).Fine(2) - y0_DTM)/dY_DTM);
    S(i).Fine(3) = DTM(i_y_fine,i_x_fine);

    A_cella = abs(dX_DTM * dY_DTM); % area della cella (= 1m x 1m = 1 m2)
    % se raggio di ricerca rdr = 5-1 = 4 m, allora la matrice sarà una 9x9
    dim_matrice = 15; % dimensione della matrice QUADRATA; DEVE AVERE DIMENSIONE DISPARI; è un parametro DA TARARE
    rdr = (dim_matrice - 1)/2; % raggio di ricerca (m), per cercare il max in una matrice 10x10 centrata nel punto considerato

    if S(i).sorgenteOcollettore == "sorgente"
        S(i).A_monte_sorg = max(max(accumulazione(i_y_inizio- rdr :i_y_inizio+ rdr , i_x_inizio- rdr :i_x_inizio+ rdr))) * A_cella / (10^4); % (ha)
    else % si mette 0 nei collettori per poter calcolare dopo l'area totale del ramo
        S(i).A_monte_sorg = 0;
    end
end

% calcolo dell'area propria di ogni ramo

for i = 1 : length(S)
    S(i).A_propria = ((A_tot_mappa - sum([S.A_monte_sorg])) / ...
        (sum([S.L]) * sum([S.n_strahler]))) * S(i).L * S(i).n_strahler; % (ha)
end

% calcolo dell'area completa dei rami sorgenti
for i = 1 : length(S)
    if S(i).sorgenteOcollettore == "sorgente"
        S(i).A_completa_sorg = S(i).A_propria + S(i).A_monte_sorg; % (ha)
    else
        S(i).A_completa_sorg = 0;
    end
end

% determinazione dell'area a monte:

rami_da_visitare = [S.fid]';
i = 1;
n_rami_rimanenti = length(S);
while ~isempty(rami_da_visitare)

    indice_ramo = find([S.fid] == rami_da_visitare(i));

    if S(indice_ramo).sorgenteOcollettore == "sorgente"
        S(indice_ramo).A_amonte = S(indice_ramo).A_completa_sorg; % (ha)

        S(indice_ramo).area_x_fi = S(indice_ramo).A_completa_sorg *  S(indice_ramo).fi; % necessario per il calcolo fi equivalente
        S(indice_ramo).area_tot = S(indice_ramo).A_completa_sorg; % necessario per il calcolo fi equivalente

        rami_da_visitare(i) = [];
        n_rami_rimanenti = n_rami_rimanenti - 1;
    else
        convergenti_in_i = [find([S.convergein] == S(indice_ramo).fid)];

        Area_a_monte_i = 0;
        area_tot_i = 0; % necessario per il calcolo fi equivalente
        area_x_fi_i = 0; % necessario per il calcolo fi equivalente
        stop = 0; % variabile che serve per capire quando si sta considerando un ramo di cui non si conosce l'area a monte

        for j = 1:length(convergenti_in_i)

            if isempty(S(convergenti_in_i(j)).A_amonte) % se non si ha l'area a monte di ogni ramo convergente in quello considerato
                stop = 1;
                break
            else
                Area_a_monte_i = Area_a_monte_i + S(convergenti_in_i(j)).A_amonte;

                area_x_fi_i = area_x_fi_i + S(convergenti_in_i(j)).area_x_fi; % necessario per il calcolo fi equivalente
                area_tot_i = area_tot_i + S(convergenti_in_i(j)).area_tot; % necessario per il calcolo fi equivalente

            end
        end
        % una volta uscito dal for, se tutti i rami convergenti in quello
        % considerato avevano l'area a monte e numero di strahler già calcolati, allora stop = 0
        % e posso calcolare l'area monte e il numero di strahler del ramo considerato; altrimenti
        % si passa al ramo successivo incrementando k

        if stop == 0
            S(indice_ramo).A_amonte = Area_a_monte_i + S(indice_ramo).A_propria;

            S(indice_ramo).area_x_fi = S(indice_ramo).A_propria * S(indice_ramo).fi + area_x_fi_i; % necessario per il calcolo fi equivalente
            S(indice_ramo).area_tot = S(indice_ramo).A_propria + area_tot_i; % necessario per il calcolo fi equivalente

            rami_da_visitare(i) = [];
            n_rami_rimanenti = n_rami_rimanenti - 1;
        else
            i = i + 1;
        end
    end
    % se sono arrivato alla fine dei rami, allora si mette k = 1 e si
    % ricomincia da capo per calcolare le aree dei rami non calcolati precedentemente;
    % se non ci sono più rami da visitare, allora esce dal while
    if i == n_rami_rimanenti+1
        i = 1;
    end
end


%% calcolo del coefficiente di afflusso equivalente
for i=1:length(S)
    if S(i).sorgenteOcollettore == "sorgente"
        S(i).fi_eq = S(i).fi;
    else
        S(i).fi_eq = (1/S(i).area_tot) * S(i).area_x_fi;
    end
end


%% calcolo della pendenza dei rami
% calcolata considerando solo i punti di inizio e di fine e la lunghezza reale del ramo L (%):

for i = 1:length(S)
    S(i).i_f_iniziofine = ((S(i).Inizio(3) - S(i).Fine(3)) / S(i).L) * 100; % pendenza (%)
end


%% controllo con le quote del DTM se i rami sono in discesa o in salita

for i = 1 : length(S)
    S(i).diff_z = S(i).Inizio(3) - S(i).Fine(3);
    if S(i).diff_z < 0
        S(i).dif_z = 'ramo IN SALITA!';
    else
        S(i).dif_z = 'ramo in discesa';
    end
end


%% determinazione della pendenza e delle quote della rete:

for i = 1 : length(rami_terminali)
    rami_adiacenti = [];
    pendenze = []; % dichiarazione di un vettore in cui si inseriscono le pendenze dei singoli rami

    indice_ramo = find([S.fid] == rami_terminali(i));

    S(indice_ramo).quota_fineramo = 0; % assegnazione ai rami terminali la quota finale = 0 m

    S(indice_ramo).quote_scavi = 3.3; % (m) parametro di TARATURA
    S(indice_ramo).quota_inizioramo = S(indice_ramo).Inizio(3) - S(indice_ramo).quote_scavi;
    S(indice_ramo).pendenze = ((S(indice_ramo).quota_inizioramo - S(indice_ramo).quota_fineramo)/ S(indice_ramo).L) * 100;

    for j = 1 : length(S)
        if S(indice_ramo).fid == S(j).convergein % se la fid del ramo che si sta considerando è = al valore di convergein di un altro ramo
            rami_adiacenti = [rami_adiacenti, S(j).fid]; % allora si aggiunge a rami_adiacenti la fid del ramo convergente
            S(j).quota_fineramo = S(indice_ramo).quota_inizioramo; % si assegna come quota finale del ramo convergente la quota iniziale del ramo in cui converge (si sta risalendo la rete dal basso verso l'alto)
        end
    end

    while ~isempty(rami_adiacenti) % finché non sono stati controllati tutti i rami non esce dal while
        indice_ramo_adiacente = find([S.fid] == rami_adiacenti(1));

        valore = 11;
        if S(indice_ramo_adiacente).Inizio(3) > valore
            S(indice_ramo_adiacente).quote_scavi = 1;
            S(indice_ramo_adiacente).quota_inizioramo = S(indice_ramo_adiacente).Inizio(3) - S(indice_ramo_adiacente).quote_scavi;
            S(indice_ramo_adiacente).pendenze = ((S(indice_ramo_adiacente).quota_inizioramo - S(indice_ramo_adiacente).quota_fineramo)/ S(indice_ramo_adiacente).L) * 100;

        else

            if S(indice_ramo_adiacente).i_f_iniziofine < 0.20
                S(indice_ramo_adiacente).quota_inizioramo = S(indice_ramo_adiacente).quota_fineramo + S(indice_ramo_adiacente).L * 0.2/100;
                S(indice_ramo_adiacente).pendenze = 0.2;
            else
                S(indice_ramo_adiacente).quota_inizioramo = S(indice_ramo_adiacente).quota_fineramo + S(indice_ramo_adiacente).L * S(indice_ramo_adiacente).i_f_iniziofine/100;
                S(indice_ramo_adiacente).pendenze = S(indice_ramo_adiacente).i_f_iniziofine;
            end
            S(indice_ramo_adiacente).quote_scavi = S(indice_ramo_adiacente).Inizio(3) - S(indice_ramo_adiacente).quota_inizioramo;

        end


        for j = 1 : length(S)
            if S(indice_ramo_adiacente).fid == S(j).convergein % se la fid del ramo che si sta considerando è = al valore di convergein di un altro ramo
                rami_adiacenti = [rami_adiacenti, S(j).fid]; % allora si aggiunge a rami_adiacenti la fid del ramo convergente
                S(j).quota_fineramo = S(indice_ramo_adiacente).quota_inizioramo; % si assegna come quota finale del ramo convergente la quota iniziale del ramo in cui converge (si sta risalendo la rete dal basso verso l'alto)
            end
        end
        rami_adiacenti(1) = []; % si toglie da questo vettore il ramo appena considerato (altrimenti non funzionerebbe il while - andrebbe in loop)
    end

    stop = 0;
end


%% disp delle statistiche sulle profondità di scavo

numeri_tra_0_e_1 = 0;
numeri_negativi = 0;
scavi_sopra_soglia = 0;
pendenza_minore_02 = 0;
for i=1:length(S)
    % quote_scavi(t) = S(t).quote_scavi;
    if S(i).quote_scavi > 0 && S(i).quote_scavi < 1
        numeri_tra_0_e_1 = numeri_tra_0_e_1 + 1;
    end
    if S(i).quote_scavi < 0
        numeri_negativi = numeri_negativi + 1;
    end
    if S(i).quote_scavi > 2.5
        scavi_sopra_soglia = scavi_sopra_soglia + 1;
    end
end

media = mean([S.quote_scavi]);
massimo = max([S.quote_scavi]);
minimo = min([S.quote_scavi]);
disp('-----------------------------------------------------------------')
disp('PROFONDITA'' DI SCAVO:')
disp([' • la profondità di scavo media è ' num2str(media) ' m'])
disp([' • la profondità di scavo massima è ' num2str(massimo) ' m'])
disp([' • la profondità di scavo minima è ' num2str(minimo) ' m'])
disp([' • ' num2str(numeri_tra_0_e_1) ' condotte richiedono uno scavo compreso tra 0 m e 1 m'])
disp([' • ' num2str(numeri_negativi) ' condotte sono sopra terra'])
disp([' • ' num2str(scavi_sopra_soglia) ' condotte richiedono uno scavo superiore a 2.5 m'])


%% diagramma a barre per visualizzare le classi di profondità di scavo

profondita_scavo = [S.quote_scavi];
edges = 0 : 0.5 : 4.5; % definizione degli intervalli di profondità scavo (bin edges)
[N, edges] = histcounts(profondita_scavo, edges); % si utilizza di histcounts per contare il numero di occorrenze in ciascun intervallo

numero_figure = numero_figure + 1;
figure(numero_figure)
bar(edges(1:end-1), N, 'histc'); % si utilizza l'inizio di ciascun bin per l'asse x
xlabel('Classi di profondità di scavo (m)'); ylabel('Numero di scavi');



ax = gca;
ax.XGrid = 'off'; % disattivazione delle linee di griglia verticali
ax.YGrid = 'on'; % attivazione le linee di griglia orizzontali

% se si utilizza il comando 'histcounts', il comportamento predefinito è di
% includere il valore di confine superiore nel bin successivo (quindi, se
% si ha un valore pari a 2.5, esso sarà incluso nel bin 2.5-3 e non in 2-2.5)

% per verificare i valori del diagramma 'bar'
contatore_profonditadiscavo = 0;
for i = 1 : length(S)
    if S(i).quote_scavi > 1.99 && S(i).quote_scavi <2.5
        contatore_profonditadiscavo = contatore_profonditadiscavo + 1;
    end
end


%% disp delle statistiche sulle pendenze delle condotte

disp('-----------------------------------------------------------------')
disp('PENDENZE DELLE CONDOTTE:')

pendenze_values = [S.pendenze];
valori_fid = [S.fid];
[~, indice_max] = max(pendenze_values); % si trova l'indice del valore massimo nella colonna L
max_fid = valori_fid(indice_max); % si trova il valore corrispondente nella colonna fid
max_pendenza = pendenze_values(indice_max); % si trova il valore massimo nella colonna L
disp([' • il ramo ' num2str(max_fid) ' è il ramo più pendente della rete, con una pendenza if = ' num2str(max_pendenza,'%.2f') ' %']);
pendenza_media = mean(pendenze_values);
disp([' • la pendenza media dei rami della rete è pari a ' num2str(pendenza_media,'%.2f') ' %'])
minimo = min([S.pendenze]);
disp([' • la pendenza minima della rete è pari a ' num2str(minimo) ' %'])


%% CALCOLI

for i = 1 : length(S)
    S(i).D = 0;
end

rami_da_visitare = [S.fid]';
i = 1;
n_rami_rimanenti = length(S);
while ~isempty(rami_da_visitare)

    indice_ramo = find([S.fid] == rami_da_visitare(i));
    i_D = 1;
    if S(indice_ramo).sorgenteOcollettore == "sorgente"

        [diametri(i_D), t, Wi, WM, w, u, QM, Q0] = calcolo_portate(S, diametri(i_D), w0, n, a, ks, indice_ramo);

        % se QM (portata idrologica) > Q0 (portata a speco pieno)
        while QM > Q0  % se QM > Q0 --> il D scelto prima non basta
            i_D = i_D + 1;
            [diametri(i_D), t, Wi, WM, w, u, QM, Q0] = calcolo_portate(S, diametri(i_D), w0, n, a, ks, indice_ramo);
        end

        % metodo iterativo per trovare teta_s = teta_d
        t0 = 4.43; % valore iniziale per ti [rad]
        deltaQ = 1.0000000; % inizializzazione di deltaQ
        t_avg = 0;
        cont = 0;
        while deltaQ >=0.003 % deltaQ < 3 l/s
            cont = cont + 1;
            if cont == 100
                disp(['l''iterazione per il ramo ' num2str(indice_ramo) ' NON ha raggiunto la convergenza!'])
                break
            end
            equation = @(N) (((N - sin(N)).^(5/3)) ./ (2 * pi * QM  / Q0 )).^(3/2) - N; % definizione della funzione che rappresenta l'equazione
            [t_avg, R, QP, Wi, WM, w, u, QM, deltaQ, t0] = iterazione(equation, t_avg, t0, diametri(i_D), ks, S, w0, n, a, indice_ramo);
        end

        gdr = 0.5 * (1 - cos(t_avg /2)); % [-] % appena deltaQ <= 3 l/s
        U = ks * sqrt(S(indice_ramo).pendenze/100) * R^(2/3); % (m/s)

        S(indice_ramo).D = diametri(i_D);
        S(indice_ramo).gdr = gdr;
        S(indice_ramo).Q0 = Q0;
        S(indice_ramo).QP = QP;
        S(indice_ramo).QM = QM;
        S(indice_ramo).U = U;
        S(indice_ramo).u = u;
        S(indice_ramo).w = w;
        S(indice_ramo).Wi = Wi;
        S(indice_ramo).WM = WM;
        S(indice_ramo).t_avg = t_avg;
        S(indice_ramo).deltaQ = deltaQ;

        rami_da_visitare(i) = [];
        n_rami_rimanenti = n_rami_rimanenti - 1;
    else

        convergenti_in_i = [find([S.convergein] == S(indice_ramo).fid)];
        stop = 0;

        for j = 1:length(convergenti_in_i)
            if S(convergenti_in_i(j)).D == 0 % se non si ha l'area a monte di ogni ramo convergente in quello considerato
                stop = 1;
                break
            end
        end

        % una volta uscito dal for, se tutti i rami convergenti in quello
        % considerato avevano l'area a monte e numero di strahler già calcolati, allora stop = 0
        % e si possono calcolare l'area monte e il numero di strahler del ramo considerato; altrimenti
        % si passa al ramo successivo incrementando k

        if stop == 0
            S(indice_ramo).D = max([S(convergenti_in_i).D]);

            % S(indice_ramo).D = max([S(convergenti_in_i).D])+ 0.1; % D_collettore > max{rami convergenti nel collettore}
            % D = S(indice_ramo).D;
            i_D = find(diametri == S(indice_ramo).D);

            % qui in n_convergenti si prendono i rami interessati di tutta la rete, non solo i convergenti al ramo i-esimo
            [n_convergenti] = trova_convergenti(indice_ramo,S);

            [diametri(i_D), t, Wi, WM, w, u, QM, Q0] = calcolo_portate(S, diametri(i_D), w0, n, a, ks, indice_ramo, n_convergenti);

            % se QM > Q0
            while QM >Q0  % se QM > Q0 il D scelto prima non basta
                i_D = i_D + 1;
                [diametri(i_D), t, Wi, WM, w, u, QM, Q0] = calcolo_portate(S, diametri(i_D), w0, n, a, ks, indice_ramo, n_convergenti);
            end

            % metodo iterativo per trovare teta_s = teta_d
            t0 = 4.43;
            deltaQ =1.0000000; % inizializzazione di deltaQ
            t_avg = 0;
            cont = 0;
            while deltaQ >=0.003 % deltaQ < 3 l/s
                cont = cont + 1;
                if cont == 100
                    disp(['l''iterazione per il ramo ' num2str(indice_ramo) ' non ha raggiunto la convergenza!'])
                    break
                end
                equation = @(N) (((N - sin(N)).^(5/3)) ./ (2 * pi * QM  / Q0 )).^(3/2) - N; % definizione della funzione che rappresenta l'equazione

                [t_avg, R, QP, Wi, WM, w, u, QM, deltaQ, t0] = iterazione(equation, t_avg, t0, diametri(i_D), ks, S, w0, n, a, indice_ramo, n_convergenti);
            end

            gdr  = 0.5 * (1 - cos(t_avg /2)); % appena deltaQ <= 3 l/s
            U = ks * sqrt(S(indice_ramo).pendenze/100) * R^(2/3);

            S(indice_ramo).D = diametri(i_D);
            S(indice_ramo).gdr = gdr;
            S(indice_ramo).Q0 = Q0;
            S(indice_ramo).QP = QP;
            S(indice_ramo).QM = QM;
            S(indice_ramo).U = U;
            S(indice_ramo).u = u;
            S(indice_ramo).w = w;
            S(indice_ramo).Wi = Wi;
            S(indice_ramo).WM = WM;
            S(indice_ramo).t_avg = t_avg;
            S(indice_ramo).deltaQ = deltaQ;

            rami_da_visitare(i) = [];
            n_rami_rimanenti = n_rami_rimanenti - 1;
        else
            i = i + 1;
        end
    end
    % se si è arrivati alla fine dei rami, allora si mette k = 1 e ricomincio
    % da capo per calcolare le aree dei rami non calcolati precedentemente;
    % se non ci sono più rami da visitare, allora esce dal while
    if i == n_rami_rimanenti+1
        i = 1;
    end
end


%% controlli su grado di riempimento gdr e velocità U:

for i = 1 : length(S)
    % grado di riempimento (gdr): 0.6 < gdr < 0.9
    if S(i).gdr >0.5 && S(i).gdr < 5
        S(i).controllo_gdr = 'OK';
    else
        S(i).controllo_gdr = 'NO!';
    end

    % velocità U: 0.5 < U < 5 m/s
    if S(i).U > 0.6 && S(i).U < 0.9
        S(i).controllo_U = '✓';
    else
        S(i).controllo_U = "NOO!";
    end
end


%% disp dei km di condotte necessari per il progetto (per ogni diametro commerciale)

disp('-----------------------------------------------------------------')
somma_lunghezze = [];
for i = 1:length(diametri)
    indici_con_D_iesimo = find([S.D] == diametri(i));
    somma_lunghezze(i) = sum([S(indici_con_D_iesimo).L]);
end

disp('LUNGHEZZE DI CONDOTTE NECESSARIE PER OGNI DIAMETRO')
larghezza = max(max(length(num2str(max(diametri), '%.1f')), 4), max(length(num2str(max(somma_lunghezze))), 3));
fprintf(['Diametro (m): ', repmat(['%', num2str(larghezza), '.2f '], 1, numel(diametri)), '\n'], diametri);
fprintf(['L_totale (m): ', repmat(['%', num2str(larghezza), '.0f '], 1, numel(somma_lunghezze)), '\n\n'], somma_lunghezze);

somma_totale_L = sum([S.L]);
disp(['L tot = ' num2str(somma_totale_L) ' m'])
somma_elle = sum(somma_lunghezze);
disp(['somma L = ' num2str(somma_elle) ' m'])
disp('-----------------------------------------------------------------')


%% diagramma a barre diametri commerciali - km di condotte necessari per il progetto
numero_figure = numero_figure + 1;
figure(numero_figure)
bar(diametri, somma_lunghezze/1000)

xlabel('Diametro delle condotte (m)'); ylabel('Lunghezza totale necessaria per il progetto (km)');
% title('Titolo');
ax = gca;
ax.XGrid = 'off'; % disattivazione delle linee di griglia verticali
ax.YGrid = 'on'; % attivazione le linee di griglia orizzontali
set(gca, 'TickLength', [0 0]);  % si imposta la lunghezza dei ticks = 0
set(gca, 'YTick', 0 : 2 : 22);
ylim([0 22])


%% plot del DTM (2 modi: 2D e 3D)

soglia = 0.50; % soglia sotto cui colorare tutto di azzurro mare (uniforme)
n_colormap = 256; % numero totale di colori nella colormap

% colormap uniforme per valori sotto la soglia (azzurro mare)
n_sotto_soglia = round(n_colormap * (soglia - min(DTM(:))) / (max(DTM(:)) - min(DTM(:))));
azzurro_mare = [0, 0.5, 1];
colormap_sotto_soglia = repmat(azzurro_mare, n_sotto_soglia, 1); % solo azzurro (mare)

% colormap gradiente per valori sopra la soglia
n_sopra_soglia = n_colormap - n_sotto_soglia;

% definizione dei colori del gradiente
verde = [0, 0.5, 0];  % verde scuro
giallo = [1, 1, 0]; % giallo
marrone = [0.6, 0.3, 0]; % marrone

% creazione del gradiente verde → giallo → marrone
% da verde a giallo
gradiente_verde_giallo = [linspace(verde(1), giallo(1), round(n_sopra_soglia/2))', ...
    linspace(verde(2), giallo(2), round(n_sopra_soglia/2))', ...
    linspace(verde(3), giallo(3), round(n_sopra_soglia/2))'];
% da giallo a marrone
gradiente_giallo_marrone = [linspace(giallo(1), marrone(1), round(n_sopra_soglia/2))', ...
    linspace(giallo(2), marrone(2), round(n_sopra_soglia/2))', ...
    linspace(giallo(3), marrone(3), round(n_sopra_soglia/2))'];

% unione dei due gradienti
colormap_gradiente = [gradiente_verde_giallo; gradiente_giallo_marrone];

% combinazione delle due colormap
colormap_comb = [colormap_sotto_soglia; colormap_gradiente];

% primo plot: mesh
numero_figure = numero_figure + 1;
figure(numero_figure)
ax1 = axes('Position', [0.1, 0.1, 0.75, 0.8]); % da decidere
mapshow(ax1, DTM, info_DTM, 'DisplayType', 'mesh');
colormap(colormap_comb);

% colorbar
cbar1 = colorbar('Location', 'eastoutside');
cbar1.Label.String = 'Valori Elevazione';

% definizione dei limiti della colorbar
cbar1.Limits = [min(DTM(:)), max(DTM(:))];

% calcolo dei ticks della colorbar (si parte da 25, incrementando di 25)
min_val = 25; % valore minimo dei tick
max_val = ceil(max(DTM(:))); % valore massimo dei tick (arrotondato)
ticks = min_val : 25 : max_val; % ticks della colorbar

% si aggiunge la soglia ai ticks se non è già inclusa
if soglia < min(ticks) || soglia > max(ticks)
    ticks = [min(ticks) - 25, ticks, max(ticks) + 25];
end

% impostazione dei ticks e delle etichette della colorbar
cbar1.Ticks = ticks;
tick_labels = arrayfun(@(x) sprintf('%d', x), cbar1.Ticks, 'UniformOutput', false);
% si aggiunge l'etichetta specifica per la soglia di 0.50 metri
if any(ticks == floor(soglia))
    tick_labels{ticks == floor(soglia)} = sprintf('%.2f', soglia);
end
cbar1.TickLabels = tick_labels;

hold on
plot([S.X], [S.Y], 'ro');

title('Rappresentazione del DTM e della rete di drenaggio urbano');
xlabel('X (m)');
ylabel('Y (m)');
xlim([info_DTM.XWorldLimits(1), info_DTM.XWorldLimits(2)]); % limiti asse x
ylim([info_DTM.YWorldLimits(1), info_DTM.YWorldLimits(2)]); % limiti asse y

box on

print('2D_DTM_e_rete', '-dpdf', '-r300');


% secondo plot: surface
numero_figure = numero_figure + 1;
figure(numero_figure)
ax2 = axes('Position', [0.1, 0.1, 0.65, 0.8]); % da decidere
mapshow(ax2, DTM, info_DTM, 'DisplayType', 'surface');
colormap(colormap_comb); % applicazione della colormap combinata
view(3);
set(ax2, 'Color', 'none'); % rimozione del piano bianco posteriore

% si aggiunge e si personalizza la colorbar
cbar2 = colorbar('Location', 'eastoutside');
% si ottiene la posizione degli assi per allineare correttamente la colorbar
pos_ax2 = get(ax2, 'Position');
cbar2.Position = [pos_ax2(1) + pos_ax2(3) + 0.04, pos_ax2(2) + 0.15, 0.02, pos_ax2(4) - 0.3];
cbar2.Label.String = 'Valori di elevazione (m)'; % etichetta della colorbar

% definizione dei limiti e dei ticks della colorbar
cbar2.Limits = [min(DTM(:)), max(DTM(:))]; % limiti della colorbar

% calcolo dei ticks della colorbar (si parte da 25, incrementando di 25)
min_val = 25; % valore minimo dei tick
max_val = ceil(max(DTM(:))); % valore massimo dei tick (arrotondato)
ticks = min_val:25:max_val; % ticks della colorbar

% si aggiunge la soglia ai ticks se non è già inclusa
if soglia < min(ticks) || soglia > max(ticks)
    ticks = [min(ticks) - 25, ticks, max(ticks) + 25];
end

% si impostano i ticks e le etichette della colorbar
cbar2.Ticks = ticks; % imposto i ticks
tick_labels = arrayfun(@(x) sprintf('%d', x), cbar2.Ticks, 'UniformOutput', false);
% si aggiunge l'etichetta specifica per la soglia
if any(ticks == floor(soglia))
    tick_labels{ticks == floor(soglia)} = sprintf('%.2f', soglia);
end
cbar2.TickLabels = tick_labels;

% title('Rappresentazione 3D del DTM');

xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
xlim([info_DTM.XWorldLimits(1), info_DTM.XWorldLimits(2)]); % limiti asse x
ylim([info_DTM.YWorldLimits(1), info_DTM.YWorldLimits(2)]); % limiti asse y

hold off


%% calcolo e disp di max, min e media delle grandezze D, gdr, U

colonne = {'D', 'gdr', 'U'};  % cell di stringhe
fprintf('Campo\t\tMinimo\t\tMassimo\t\tMedio\n');
fprintf('_________________________________________________________\n');
field_data = [];
for i = 1:length(colonne)
    field_data = [S.(colonne{i})]; % si ottengono i dati riferiti al campo corrente
    min_val = min(field_data);
    max_val = max(field_data);
    mean_val = mean(field_data);
    fprintf('%s\t\t%.2f\t\t%.2f\t\t%.2f\n', colonne{i}, min_val, max_val, mean_val);
end

numero_figure = numero_figure + 1;
figure(numero_figure)


%% plot della rete + numero del ramo + colorazione dei rami a seconda di un parametro (valori discreti - legenda semplice)
colonne = {'sorgenteOcollettore', 'id_sottobacino', 'n_strahler', 'D'};
for i = 1 : length(colonne)
    plot_rete_valori_discreti(colonne{i}, numero_figure, S, orto_merged, R_merged)
    numero_figure = numero_figure + 1;
end



%% plot dlla rete + scrittura di un determinato valore sui rami + colorazione dei rami a seconda di quel valore (valori continui - colorbar)
colonne = {'L', 'fi', 'quote_scavi', 'pendenze', 'D', 'gdr', 'U'};
for i = 1 : length(colonne)
    plot_rete_valori_continui(colonne{i}, numero_figure, S, orto_merged, R_merged)
    numero_figure = numero_figure + 1;
end

%% riorganizzazione dell'ordine delle colonne nella struct S
% S = orderfields(S, {'fid','id_sottobacino','sorgenteOcollettore',...
%     'convergein','X', 'Y', 'Z', 'Inizio', 'Fine', 'L', 'diff_z', 'dif_z', 'i_f_iniziofine',...
%     'n_strahler', 'fi', 'A_monte_sorg', 'A_propria', 'A_completa_sorg', 'area_tot', 'A_amonte',...
%     'area_x_fi', 'fi_eq', 'quota_inizioramo', 'quota_fineramo', 'quote_scavi',...
%     'pendenze', 'controllo_singoli_tratti', 'D', 'gdr',...
%     'Q0', 'QP', 'QM', 'U', 'u', 'w', 'Wi', 'WM', 't_avg', ...
%     'deltaQ', 'controllo_gdr', 'controllo_U'});


%% si cancella dal workspace tutto ciò che non serve più

% clear A_monte_sorg Area_a_monte_k campi_da_estrarre contatore dim_matrice i i_x_fine i_x_inizio...
%     i_y_fine i_y_inizio ix iy j k L_tmp max_difflungh rami_adiacenti_i rami_visitati rdr tmp


%% funzioni

% • calcolo_portate
% • trova_convergenti
% • iterazione: la funzione 'iterazione' permette di trovare deltaQ. le
%               iterazioni si fermano non appena il valore di detlaQ < 3 l/s
% • plot_rete_valori_continui
% • plot_rete_valori_discreti
% • accessori_plot


%% funzione: calcolo_portate

function [D, t, Wi, WM, w, u, QM, Q0] = calcolo_portate(S, D, w0, n, a, ks, i, n_convergenti)
if S(i).sorgenteOcollettore == "sorgente"
    gdr = 0.8; % h/D: grado di riempimento (ottimale) (-)
    t = 2 * acos(1 - 2*gdr); % theta - θ (rad)
    omega = ((D ^2)/8)*(t - sin(t)); % (m2)
    Wi  = S(i).L * omega; % (m3)
    WM  = w0 * S(i).A_completa_sorg + Wi ; % (m3)
    w  = double(WM / (S(i).A_completa_sorg*(10^4))); % volume di invaso specifico (m3/m2)
    u = 2168 * (n*(S(i).fi * a * (10^(-3)))^(1/n)) / (w^((1/n)-1)); % (l/(s ha))
    QM  = double((u  * S(i).A_completa_sorg)/(1000)); % portata idrologica (m3/s)
    Q0 = double((ks * sqrt(S(i).pendenze/100)*((D/4)^(2/3))*(pi*D^2)/4)); % portata idraulica - speco pieno (m3/s)
else
    % calcolo di QM e Q0
    gdr = 0.8;
    S(i).gdr = gdr;
    t = 2 * acos(1 - 2*gdr);
    omega = ((D ^2)/8)*(t - sin(t));
    Wi  = S(i).L * omega;
    S(i).Wi = Wi;

    WM = 0;
    n_convergenti = find(ismember([S.fid], n_convergenti));
    for j = 1 : length(n_convergenti)
        WM = WM + S(n_convergenti(j)).Wi;
    end
    WM  = WM + Wi + w0 * S(i).A_amonte;
    w  = double(WM  / (S(i).A_amonte*(10^4)));
    u = 2168 * (n * (S(i).fi_eq * a * (10^(-3)))^(1/n)) / (w^((1/n) - 1));
    QM  = double((u  * S(i).A_amonte)/(1000));
    Q0 = double((ks * sqrt(S(i).pendenze/100)*((D/4)^(2/3))*(pi*D^2)/4));
end
end


%% funzione: trova_convergenti

function [n_convergenti] = trova_convergenti(i,S)
n_convergenti = [];
vettore_appoggio = [];
% inizializzazione del vettore_appoggio con il ramo iniziale
vettore_appoggio = [vettore_appoggio; S(i).fid];
while ~isempty(vettore_appoggio)
    % si estrae il prossimo ramo da elaborare
    ramo_corrente= vettore_appoggio(end);
    vettore_appoggio = vettore_appoggio(1:end-1);  % rimozione dell'elemento superiore dalla pila
    % si trovano tutti i rami che convergono al ramo i-esimo

    indice_ramo = find([S.fid] == ramo_corrente);

    rami_conv_i = [];
    for j = 1 : length(S)
        % si trovano tutti i rami convergenti in quello che si sta considerando:
        if S(indice_ramo).fid == S(j).convergein
            rami_conv_i = [rami_conv_i; S(j).fid];
        end
    end

    % si aggiungono i rami convergenti al ramo i-esimo (rami_conv_i) a n_convergenti
    n_convergenti = [n_convergenti; rami_conv_i];
    % si aggiungono i rami convergenti al ramo i-esimo al vettore di appoggio
    vettore_appoggio = [vettore_appoggio; rami_conv_i];
end
end


%% funzione: iterazione

function [t_avg, R, QP, Wi, WM, w, u, QM, deltaQ,t0] = iterazione(equation, t_avg, t0, D, ks, S, w0, n, a, i, n_convergenti)
if S(i).sorgenteOcollettore == "sorgente"
    options = optimoptions('fsolve','Display','off');
    ti_solution  = fsolve(equation, t0, options); % risoluzione dell'equazione con fsolve
    t_avg  = mean([t0 ti_solution ]); % θ medio (rad)
    t0 = t_avg ;
    omega  = ((D ^2)/8)*(t_avg  - sin(t_avg )); % (m2)
    B  = t_avg  * D /2; % (m)
    R  = omega  / B ; % (m)
    QP  = ks * sqrt(S(i).pendenze/100) * omega  * R ^(2/3); % (m3/s)
    Wi  = S(i).L * omega ; % (m3)
    WM  = w0 * S(i).A_completa_sorg + Wi ; % (m3)
    w  = double(WM  / (S(i).A_completa_sorg * (10^4))); % (m3/m2)
    u = 2168 * (n*(S(i).fi * a * (10^(-3)))^(1/n)) / (w^((1/n)-1)); % (l/(s ha))
    QM  = double((u  * S(i).A_completa_sorg)/1000); % (m3/s)
    deltaQ  = abs(QP -QM ); % = abs(QM - QP) (m3/s)

else
    options = optimoptions('fsolve','Display','off');
    ti_solution  = fsolve(equation, t0, options); % risoluzione dell'equazione con fsolve
    t_avg  = mean([t0 ti_solution ]);
    t0 = t_avg;
    omega  = ((D ^2)/8) * (t_avg - sin(t_avg));
    B = t_avg * D/2;
    R = omega / B;
    QP = ks * sqrt(S(i).pendenze/100) * omega  * R ^(2/3);
    Wi = S(i).L * omega;

    WM = 0;
    n_convergenti = find(ismember([S.fid], n_convergenti));
    for j = 1 : length(n_convergenti)
        WM = WM + S(n_convergenti(j)).Wi;
    end
    WM  = WM + Wi + w0 * S(i).A_amonte;
    w  = WM  / (S(i).A_amonte * (10^4));
    u = 2168 * (n * (S(i).fi_eq * a * (10^(-3)))^(1/n)) / (w^((1/n)-1));
    QM  = double((u * S(i).A_amonte)/(1000));
    deltaQ  = abs(QP -QM );

end
end


%% funzione: plot della rete + colorazione dei rami a seconda di valori continui

function plot_rete_valori_continui(colonna, numero_figure, S, orto_merged, R_merged)

figure(numero_figure)
hold on;

% mapshow della rete con ortofoto in trasparenza (= 0.40)
mapshow(orto_merged, R_merged, 'DisplayType', 'image', 'AlphaData', 0.3);

colori = turbo(64);  % definizione di una colormap di 64 colori

% si impostano i limiti per la normalizzazione in base alla colonna
if strcmp(colonna, 'gdr') || strcmp(colonna, 'fi')
    % si forzano i limiti della colorbar tra 0 e 1
    min_valore_colonna = 0;
    max_valore_colonna = 1;
    valori_colonna = [S.(colonna)];
    valori_colonna_normalizzati = valori_colonna;  % i valori sono già tra 0 e 1
else
    % limiti dinamici (da min a max) per le altre colonne
    min_valore_colonna = min([S.(colonna)]);
    max_valore_colonna = max([S.(colonna)]);
    valori_colonna = [S.(colonna)];
    valori_colonna_normalizzati = (valori_colonna - min_valore_colonna) / (max_valore_colonna - min_valore_colonna);  % normalizzazione tra 0 e 1
end

% ciclo sulla rete e plot con il colore corrispondente
for i = 1:length(S)
    x = S(i).X; y = S(i).Y;
    indice_colore = floor(valori_colonna_normalizzati(i) * (size(colori, 1) - 1)) + 1;
    % ci si assicura che l'indice del colore sia all'interno dei limiti
    indice_colore = max(1, min(indice_colore, size(colori, 1)));
    % si ottiene il colore corrispondente dalla colormap
    colore_corrente = colori(indice_colore, :);
    plot(x, y, 'Color', colore_corrente, 'LineWidth', 1.25);
    hold on;
end

% si imposta la colormap e si aggiunge la colorbar
colormap(turbo); % si usa la colormap predefinita 'turbo'
c = colorbar;
c.Label.String = ['Valore di ' regexprep(colonna,'_',' ')];
% si impostano i limiti della colorbar
if strcmp(colonna, 'gdr') || strcmp(colonna, 'fi')
    % se la colonna è 'gdr' o 'fi', si impostano i limiti della colorbar tra 0 e 1
    caxis([0, 1]);
else
    % altrimenti si impostano i limiti della colorbar dinamicamente (dal min al max)
    caxis([min_valore_colonna, max_valore_colonna]);
end

title(['Rete colorata secondo il valore di ' regexprep(colonna,'_',' ')]);

accessori_plot(S);
end


%% funzione per plottare la rete e colorarla a seconda di valori discreti

function plot_rete_valori_discreti(colonna, numero_figure, S, orto_merged, R_merged)

figure(numero_figure)
hold on;

% mapshow della rete con ortofoto in trasparenza (= 0.40)
mapshow(orto_merged, R_merged, 'DisplayType', 'image', 'AlphaData', 0.3);

% plot della rete sopra l'ortofoto
valori_colonna = [S.(colonna)];
valori_unici = unique(valori_colonna);
colori = lines(length(valori_unici));

% si crea la legenda solo per la colonna 'sorgenteOcollettore'
if strcmp(colonna, 'sorgenteOcollettore')
    % si ordinano sorgente, collettore, terminale
    ordine_predefinito = {'sorgente', 'collettore', 'terminale'};
    [~, idx_ordinato] = ismember(ordine_predefinito, valori_unici);
    valori_unici = valori_unici(idx_ordinato);

    % si ordinano anche i colori seguendo lo stesso ordine
    colori = colori(idx_ordinato, :);

    % si invertono i colori di 'sorgente' e 'collettore'
    idx_sorgente = find(strcmp(ordine_predefinito, 'sorgente'));
    idx_collettore = find(strcmp(ordine_predefinito, 'collettore'));
    temp = colori(idx_sorgente, :);
    colori(idx_sorgente, :) = colori(idx_collettore, :);
    colori(idx_collettore, :) = temp;

end

for i = 1:length(S)
    colore_corrente = colori(valori_unici == valori_colonna(i), :);
    plot(S(i).X, S(i).Y, 'Color', colore_corrente, 'LineWidth', 1.25);

    hold on;
end

title(['Rete colorata secondo il valore di ' regexprep(colonna,'_',' ')]);

% si crea la legenda
legend_str = cell(1, length(valori_unici));
for i = 1:length(valori_unici)
    legend_str{i} = [num2str(valori_unici(i))];
end

linee_legenda = zeros(1, numel(valori_unici));
for i = 1 : numel(valori_unici)
    linee_legenda(i) = plot([NaN, NaN], 'Color', colori(i,:));
    set(linee_legenda(i), 'LineWidth', 2);
    hold on;
end
legend(linee_legenda, legend_str, 'Location', 'northeast');
legend('boxoff');

accessori_plot(S)

end


%% funzione: accessori_plot

function accessori_plot(S)

axis equal
box on

% calcolo dei limiti degli assi
x_min_rete = min([S.X]);
x_max_rete = max([S.X]);
y_min_rete = min([S.Y]);
y_max_rete = max([S.Y]);

% calcolo dei margini (si possono modificare i valori)
margine_x = 0.15 * (x_max_rete - x_min_rete);
margine_y = 0.05 * (y_max_rete - y_min_rete);

% impostazione dei nuovi limiti degli assi con i margini
xlim([x_min_rete - margine_x, x_max_rete + margine_x]);
ylim([y_min_rete - margine_y, y_max_rete + margine_y]);

xlabel('X (m)');
ylabel('Y (m)');

end

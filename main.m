%Preuzimanje podataka sa linka
data = [[8.65124763, 15.05801139]; [.12218972, 18.71889232]; 
[.7461504e-1, 19.08078948]; [-1.559782888, 19.98371172]; 
[-9.827561842, 23.06848486]; [.32311207, 19.14564218]; 
[-2.391913592, 20.30217888]; [3.94617018, 17.43607480]; 
[-2.830376788, 19.63133334]; [-1.358328420, 20.05752868]; 
[6.23095437, 16.43323636]; [5.61148003, 17.81101690]; 
[-2.625640890, 20.88674163]; [-2.690509364, 20.36160632]; 
[-8.857403420, 23.84930474]; [-2.780958614, 19.81177275]; 
[.49220217, 18.73849995]; [-6.543130054, 21.51666222]; 
[6.08427501, 16.97527217]; [-6.242938578, 21.72698882]; 
[-6.514007062, 21.45289388]; [7.76133332, 16.75089278]; 
[-9.054911088, 22.54190263]; [5.03304675, 16.82722949]; 
[-9.869721387, 23.70277711]; [-7.770840930, 22.16242030]; 
[-1.948518384,20.06695018]; [4.00205027, 17.47005315]; 
[5.17321040, 17.55005418];
[-2.076084772, 19.72089223]; [-1.550062590, 20.21390427]; 
[-6.143925098, 21.65680593]; [5.60762551, 16.84921688]; 
[6.64081052, 17.04952989]; 
[2.71911796, 17.75913317]; [.49583125, 18.89613761]; 
[5.12446722, 16.91928245]; [-1.852442276, 20.25908204]; 
[6.90280682, 16.84990814]; [-2.532240154, 20.72224765]; 
[1.32201814, 18.91342503]; [.93254506, 18.90640481]; 
[-8.184966738, 29.44895089]; [-9.307860036, 27.09063691];
[-9.710001259, 32.12397498]; [-6.499868778, 28.48508601]; 
[7.492109388, 11.06831268]; [5.295254201, 9.403916851];
[8.105629482, 7.166928438]; [9, 215.6599449]; 
[-10, -191]; [-8, 200]; [9, 220]; [1, -100]];

x = data(:,1);
y = data(:,2);
n = length(x);

%Priprema podataka za rad
m = 1 %Stepen aproksimacionog polinoma
m = m+1;
p = ones(n,1);

%Nameštanje slike
figure;
screen_size = get(0, 'ScreenSize');
set(gcf, 'Position', [(screen_size(3)-1080)/2, (screen_size(4)-540)/2, 1080, 540]);

%Pre RLS
c0 = leastSquares(x, y, m, p);

disp(strjoin(["Aproksimacioni polinom pre RLS je:",  polyString(c0)], " "));
plotting(x,y,c0,1);

%Posle RLS
h = diag(x*((x'*x)\x'));

maxIter = 100;
iter = 0;
c1 = c0+ones(m,1);
while round(c1, 5) ~= round(c0, 5)
    iter = iter + 1;
    if iter > maxIter
        disp('Maximum iterations reached. Exiting...');
        break;
    end
    c1 = c0;
    p = adjustedWeights(x, y, c0, h);
    c0 = leastSquares(x, y, m, p);
end

disp(strjoin(["Aproksimacioni polinom nakon RLS je:",  polyString(c0)], " "));
plotting(x,y,c0,2);

%
%

function plotting(x, y, c, i) %F-ja za iscrtavanje grafika
    subplot(1,2,i);

    scatter(x, y, 'DisplayName', 'Подаци');
    hold on;

    xValues = linspace(-15, 15, 150);
    yValues = polyval(c, xValues);
    plot(xValues, yValues, 'DisplayName', 'Апркосимациони полином');
    
    set(gca,'position',[i*0.5-0.475, 0.05, 0.45, 0.9]);
    ylim([5, 30]); %Ako se vide outlieri onda nećemo lepo videti koliko nam je dobra aproksimacija na relevantnim vrednostima
    xlim([-11, 11]); 
   
    if i==1
        title('Метод најмањих квадрата');
    else
        title('Робустан метод најмањих квадрата');
    end
    xlabel('x');
    ylabel('y');
    legend;

    hold off;
end

function pS = polyString(c) %Postoji ugrađena f-ja poly2str, ali za nju treba poseban toolbox
    degree = length(c) - 1;
    terms = {};
    
    for i = 1:(degree + 1)
        power = degree - i + 1;
        if c(i) ~= 0
            if power == 0
                terms{end+1} = num2str(c(i));
            elseif power == 1
                terms{end+1} = [num2str(c(i)), '*x'];
            else
                terms{end+1} = [num2str(c(i)), '*x^', num2str(power)];
            end
        end
    end

    if isempty(terms)
        pS = '0';
    else
        pS = strjoin(terms, ' + ');
    end
end

function sP = scalProd(p, g, f) %F-ja koja računa naš skalarni proizvod
    sP = (p .* g)' * f;
end

function lS = leastSquares(x, y, m, p) %F-ja koja vraća koeficijente aproksimacionog polinoma preko formule sa časa
    X = ones(length(x),1);
    for i = 2 : m
        X = [X(:,i-1) .* x, X];
    end

    lS = (X'*diag(p)*X) \ (X'*diag(p)*y);
end

function aW = adjustedWeights(x, y, c, h) %Funkcija koja vraća podešene težine koristeći datu formulu
    r = (y-polyval(c, x)).^2;
    radj = r./sqrt(1-h);

    MAD = median(abs(r - median(r)));
    u = min(radj*0.6745/(4.685*MAD),1);
    
    aW = (1-u.^2).^2;
end



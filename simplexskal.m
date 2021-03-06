load(input('Ange datafil: ','s'))
tic;
[m,n] = size(A);

% Skapa nix, dvs indexvektorn för ickebasvariabler
nix = setdiff([1:n],bix);

% Skapa initial partition
B = A(:,bix);
N = A(:,nix);
cB = c(bix,:);
cN = c(nix,:);
opt=0;
iter = 0;
while opt==0
    iter = iter + 1;

    % Beräkna reducerad kostnad och högerled
    % Beräkna mest negativ reducerad kostnad, rc_min,
    % och index för inkommande variabel, inkix
    xB = B\b; % Lösa ut x ur Bx = b.
    
    y = transpose(transpose(cB) * inv(B)); % Beräkna skuggpriser A^T * y
    
    z = transpose(b) * y; % Beräkna målfunktionsvärdet
    
    cHattN = cN - transpose(N) * y; % Beräkna cHattN
    
    rc_min = cHattN(1); % Börja med att sätta rc_min till första värdet i cHattN
    
    inkix = 1; % vilket innebär att inkommande variabel för tillfället sätts till 1.
    
    for r = 1:size(cHattN) %iterera över hela cHattN
        if rc_min > cHattN(r)%Kolla om ett element i cHattN är mindre än vårt nuvarande för rc_min
            rc_min = cHattN(r); %Om det stämmer sätt det till rc_min
            inkix = r; %och sätt inkix så att den korresponderar med rätt index i cHattN
        end
    end

    
    if rc_min >= -1.0E-10
        opt=1;
        disp('Optimum');
    else
        % Beräkna inkommande kolumn, a
        a = B \ N(:,inkix);
        if max(a) <= 0 
            disp('Obegränsad lösning');
            return;
        else
            % Bestäm utgående variabel, utgix
            minValue = inf; % Börja med att sätta minsta värdet till inf
            for i = 1:size(a) % iterera över inkommande kolumn
                if xB(i) ~= 0 % kolla så att vi dividerar på noll
                    if xB(i)/a(i) > 0 % kollar att kvoten av xB(i)/a(i) inte är negativ (ej tillåten)
                        if xB(i)/a(i) < minValue % om xB(i)/a(i) är mindre än minValue vill vi:
                            minValue = xB(i)/a(i); % Ansätt minValue till det värdet 
                            utgix = i; % och spåra vilken variabel som är utgående genom att sätta utgix till i
                        end
                    end
                end
            end
                        
            
            %ratios = xB./a;
            %ratios(ratios < 0) = inf;
            %[utg, utgix] = min(ratios);


            fprintf('Iter: %d, z: %f, rc_min: %f, ink: %d, utg: %d \n',iter,z,rc_min,nix(inkix),bix(utgix));

            % Konstruera ny partitionering mha ink och utg

            temporary = bix(utgix);
            bix(utgix) = nix(inkix);
            nix(inkix) = temporary;

            % Skapa initial partition
            B = A(:,bix);
            N = A(:,nix);
            cB = c(bix,:);
            cN = c(nix,:);
        end
    end
end

toc
fprintf('z: %f\n',z);
x = zeros(n,1);
x(bix) = xB;
fprintf('sum(x-xcheat): %f\n',sum(x-xcheat));
fprintf('z-zcheat: %f\n',z-zcheat);


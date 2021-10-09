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
xN = nix;
opt=0;
iter = 0;
while opt==0
    iter = iter + 1;

    % Beräkna reducerad kostnad och högerled
    % Beräkna mest negativ reducerad kostnad, rc_min,
    % och index för inkommande variabel, inkix
    xB = inv(B) * b;
    xN = 0;
    y = transpose(transpose(cB) * inv(B));
    z = transpose(b) * y;
    
    cHatN = cN - transpose(N) * y;
    [rc_min, inkix] = min(cHatN);
    
    %for r = 1:size(cN)
    %    if r == 1
    %        rc_min = cN(r);
    %        inkix = r;
    %    end
    %    if rc_min > cN(r)
    %        rc_min = cN(r);
    %        inkix = r;
    %    end
    %end
    if rc_min >= -1.0E-10
        opt=1;
        disp('Optimum');
    else
        % Beräkna inkommande kolumn, a
        % --------
        a = N(:,inkix);
        %a = B \ N(:,inkix);
        if max(a) <= 0 
            disp('Obegränsad lösning');
            return;
        else
            % Bestäm utgående variabel, utgix
            % -------
            %for r = 1:size(a)
                %if a(r) == 0
                %    continue;
                %end
            %    if r == 1
            %        quota = b(r)/a(r);
            %        utgix = r;
            %    end
            %    if quota > b(r)/a(r)
            %        quota = b(r)/a(r);
            %        utgix = r;
            %    end
            %end
            
            ratios = xB./a;
            ratios(ratios < 0) = inf;
            [utg, utgix] = min(ratios);


            fprintf('Iter: %d, z: %f, rc_min: %f, ink: %d, utg: %d \n',iter,z,rc_min,nix(inkix),bix(utgix));

            % Konstruera ny partitionering mha ink och utg
            % --------
            temporary = bix(utgix);
            bix(utgix) = nix(inkix);
            nix(inkix) = temporary;

            xB = 0;

            temporary = B(:,utgix);
            B(:,utgix) = N(:,inkix);
            N(:,inkix) = temporary;

            temporary = cB(utgix);
            cB(utgix) = cN(inkix);
            cN(inkix) = temporary;

        end
    end
end

toc
fprintf('z: %f\n',z);
x = zeros(n,1);
x(bix) = xB;
fprintf('sum(x-xcheat): %f\n',sum(x-xcheat));
fprintf('z-zcheat: %f\n',z-zcheat);


function u = deterministicFAB_LocallyAdaptiveStepSize(u,t_max,h_1,h_2,g)

[m,n] = size(u);
U = ones(m+2,n+2);
%%Neumann boundary conditions:
U(2:m+1,2:n+1) = u;
U(1,1) = u(1,1); U(1,2:n+1) = u(1,1:n); U(1,n+2) = u(1,n);
U(2:m+1,1) = u(1:m,1); U(m+2,1) = u(1,m);
U(m+2,2:n+1) = u(m,1:n); U(m+2,n+2) = u(m,n);
U(2:m+1,n+2) = u(1:m,n);

% Initialisation:
%% Time levels for pixels:
time = zeros(m,n);

%% Diffusivity computation:
G = ones(m+2,n+2);
for i = 2:m+1
    for j = 2:n+1
        A = max(((U(i+1,j)-U(i,j))/h_1)*((U(i,j)-U(i-1,j))/h_1),0);
        B = max(((U(i,j+1)-U(i,j))/h_2)*((U(i,j)-U(i,j-1))/h_2),0);
        G(i,j) = g(A+B);
    end 
end

%% Two-pixel flow computation:
% (:,:,1) is the (i,j+1) neighbour, (:,:,2) is the (i+1,j) neighbour.
Phi = ones(m+2,n+2,2);  
for i = 2:m+1
    for j = 2:n+1
        Phi(i,j,1) = ((G(i,j) + G(i,j+1))/2)*((U(i,j+1)-U(i,j))/h_1^2);
        Phi(i,j,2) = ((G(i,j) + G(i+1,j))/2)*((U(i+1,j)-U(i,j))/h_2^2);
    end
end

%% Pixwlwise total velocity computation:
u_dt = ones(m,n);
for i = 2:m+1
    for j = 2:n+1
        u_dt(i-1,j-1) = ((G(i+1,j) + G(i,j))/2)*((U(i+1,j)-U(i,j))/(h_1^2)) - ...
               ((G(i,j) + G(i-1,j))/2)*((U(i,j)-U(i-1,j))/(h_1^2)) + ...
               ((G(i,j+1) + G(i,j))/2)*((U(i,j+1)-U(i,j))/(h_2^2)) - ...
               ((G(i,j) + G(i,j-1))/2)*((U(i,j)-U(i,j-1))/(h_2^2));
    end 
end
%for i = 1:m
%    for j = 1:n
%        u_dt(i,j) = Phi(i+1,j+1,2) - Phi(i,j+1,2) + Phi(i+1,j+1,1) - Phi(i+1,j,1);
%    end
%end

%% Flow expiry times:
T = zeros(m,n,2);
T(:,:,:) = t_max;
for i = 1:m
    for j = 1:n
        if j~=n
            i_tag = i; j_tag = j+1;
            if u(i,j)~=u(i_tag,j_tag) && sign(u_dt(i,j)-u_dt(i_tag,j_tag))==-sign(u(i,j)-u(i_tag,j_tag))...
                    && abs(u(i,j)-u(i_tag,j_tag))<T(i,j,1)*abs(u_dt(i,j)-u_dt(i_tag,j_tag))
                T(i,j,1) = -(u(i,j)-u(i_tag,j_tag))/(u_dt(i,j)-u_dt(i_tag,j_tag));
            end 
        end
        if i~=m
            i_tag = i+1; j_tag = j;
            if u(i,j)~=u(i_tag,j_tag) && sign(u_dt(i,j)-u_dt(i_tag,j_tag))==-sign(u(i,j)-u(i_tag,j_tag))...
                    && abs(u(i,j)-u(i_tag,j_tag))<T(i,j,2)*abs(u_dt(i,j)-u_dt(i_tag,j_tag))
                T(i,j,2) = -(u(i,j)-u(i_tag,j_tag))/(u_dt(i,j)-u_dt(i_tag,j_tag));
            end
        end
    end
end

%% Pixel expiry times:
T_pixel = zeros(m,n);
T_pixel(:) = t_max;
for i = 1:m
    for j = 1:n
        if u_dt(i,j)>0
            [i_tag,j_tag] = FindMaxNeighbour(u,i,j);
            if (u(i,j) + T_pixel(i,j)*u_dt(i,j))>u(i_tag,j_tag)
                T_pixel(i,j) = (u(i_tag,j_tag)-u(i,j))/u_dt(i,j);
            end
        end
        if u_dt(i,j)<0
            [i_tag,j_tag] = FindMaxNeighbour(-u,i,j);
            if (u(i,j) + T_pixel(i,j)*u_dt(i,j))<u(i_tag,j_tag)
                T_pixel(i,j) = (u(i_tag,j_tag)-u(i,j))/u_dt(i,j);
            end
        end
    end
end
        
%% Priority queue (highest priority with lowest expiry time):
Q(:,:,1) = 1.-T(:,:,1);
Q(:,:,2) = 1.-T(:,:,2);
Q(:,:,3) = 1.-T_pixel;
data(:,:,1:2) = T;
data(:,:,3) = T_pixel;

% Asynchronous update:
%timePrev = time;
cntWhile_nnz = 0; 
while nnz(time-t_max)~=0
    cntWhile_nnz = cntWhile_nnz + 1;
    cntNumPixUpdate = 0; 

    %% Priority-based selection:
    T_star = datasample(data(:),1,'Weights',Q(:));
    [i_f,j_f] = find(data==T_star);
    [size_1,~] = size(i_f);
    for i_1 = 1:size_1
        clear I_1 I_2;
        I_1(1,1) = i_f(i_1,1);
        if j_f(i_1,1)<=n
                I_1(1,2) = j_f(i_1,1); I_1(2,1) = I_1(1,1); I_1(2,2) = I_1(1,2) + 1;

                I_2(1:2,:) = I_1; I_2(3,1) = I_1(1,1); I_2(3,2) = I_1(2,2) + 1;
                I_2(4,1) = I_1(1,1); I_2(4,2) = I_1(1,2) - 1; 
                I_2(5,1) = I_1(1,1) - 1; I_2(5,2) = I_1(1,2); I_2(6,1) = I_1(1,1) + 1; I_2(6,2) = I_1(1,2);
                I_2(7,1) = I_1(1,1) - 1; I_2(7,2) = I_1(2,2); I_2(8,1) = I_1(1,1) + 1; I_2(8,2) = I_1(2,2); 
            k = 1;
        else
            if j_f(i_1,1)>n && j_f(i_1,1)<=2*n
                    I_1(1,2) = j_f(i_1,1) - n; I_1(2,1) = I_1(1,1) + 1; I_1(2,2) = I_1(1,2);

                    I_2(1:2,:) = I_1; I_2(3,1) = I_1(1,1); I_2(3,2) = I_1(1,2) + 1;
                    I_2(4,1) = I_1(1,1); I_2(4,2) = I_1(1,2) - 1;
                    I_2(5,1) = I_1(2,1); I_2(5,2) = I_1(2,2) + 1; I_2(6,1) = I_1(2,1); I_2(6,2) = I_1(2,2) - 1;
                    I_2(7,1) = I_1(1,1) - 1; I_2(7,2) = I_1(1,2); I_2(8,1) = I_1(2,1) + 1; I_2(8,2) = I_1(2,2);
            k = 2;
            else
                    I_1(1,2) = j_f(i_1,1) - 2*n; I_1(2,1) = I_1(1,1); I_1(2,2) = I_1(1,2) + 1;
                    I_1(3,1) = I_1(1,1); I_1(3,2) = I_1(1,2) - 1; I_1(4,1) = I_1(1,1) + 1;
                    I_1(4,2) = I_1(1,2); I_1(5,1) = I_1(1,1) - 1; I_1(5,2) = I_1(1,2);

                    I_2(1:5,:) = I_1; I_2(6,1) = I_1(1,1); I_2(6,2) = I_1(1,2) + 2; I_2(7,1) = I_1(1,1); I_2(7,2) = I_1(1,2) - 2;
                    I_2(8,1) = I_1(1,1) + 1; I_2(8,2) = I_1(1,2) + 1; I_2(9,1) = I_1(1,1) + 1; I_2(9,2) = I_1(1,2) - 1;
                    I_2(10,1) = I_1(1,1) + 2; I_2(10,2) = I_1(1,2); I_2(11,1) = I_1(1,1) - 1; I_2(11,2) = I_1(1,2) + 1;
                    I_2(12,1) = I_1(1,1) - 1; I_2(12,2) = I_1(1,2) - 1; I_2(13,1) = I_1(1,1) - 2; I_2(13,2) = I_1(1,2);

                k = 3;
            end
        end

        %% Pixel updates:
        [s1,~] = size(I_2);
        for i = 1:s1
            i_s = I_2(i,1); j_s = I_2(i,2);
            if i_s>0 && i_s<=m && j_s>0 && j_s<=n
                U(i_s+1,j_s+1) = U(i_s+1,j_s+1) + (T_star - time(i_s,j_s))*u_dt(i_s,j_s); 
                if (time(i_s,j_s) ~= T_star )
                    time(i_s,j_s) = T_star; 
                    cntNumPixUpdate = cntNumPixUpdate + 1; 
                end
                %u = U(2:m+1,2:n+1);
                %U(1,1) = u(1,1); U(1,2:n+1) = u(1,1:n); U(1,n+2) = u(1,n);
                %U(2:m+1,1) = u(1:m,1); U(m+2,1) = u(1,m);
                %U(m+2,2:n+1) = u(m,1:n); U(m+2,n+2) = u(m,n);
                %U(2:m+1,n+2) = u(1:m,n);
            end
        end

        %% Diffusivity computation:
        [s1,~] = size(I_1);
        for i = 1:s1
            i_s = I_1(i,1); j_s = I_1(i,2);
            if i_s>0 && i_s<=m && j_s>0 && j_s<=n
                a = max(((U(i_s+2,j_s+1) - U(i_s+1,j_s+1))*(U(i_s+1,j_s+1) - U(i_s,j_s+1)))/(h_1^2),0);
                b = max(((U(i_s+1,j_s+2) - U(i_s+1,j_s+1))*(U(i_s+1,j_s+1) - U(i_s+1,j_s)))/(h_2^2),0);
                G(i_s+1,j_s+1) = g(a+b);       
            end
        end

        %% Two-pixel flow computation:
        i_s = I_1(1,1); j_s = I_1(1,2); 
        i_ss = I_1(2,1); j_ss = I_1(2,2);
        if i_s>0 && i_s<=m && j_s>0 && j_s<=n && i_ss>0 && i_ss<=m && j_ss>0 && j_ss<=n
            if k==1                                     % (i,j+1) neighbour.
                Phi(i_s+1,j_s+1,1) = ((G(i_s+1,j_s+1) + G(i_ss+1,j_ss+1))/2)*((U(i_ss+1,j_ss+1) - U(i_s+1,j_s+1))/h_1^2);
            end
            if k==2                                     % (i+1,j) neighbour.
                Phi(i_s+1,j_s+1,2) = ((G(i_s+1,j_s+1) + G(i_ss+1,j_ss+1))/2)*((U(i_ss+1,j_ss+1) - U(i_s+1,j_s+1))/h_2^2);
            end
            if k==3
               Phi(i_s+1,j_s+1,1) = ((G(i_s+1,j_s+1) + G(i_ss+1,j_ss+1))/2)*((U(i_ss+1,j_ss+1) - U(i_s+1,j_s+1))/h_1^2);
               i_ss = I_1(3,1); j_ss = I_1(3,2);        % (i,j-1) neighbour.
               if i_ss>0 && i_ss<=m && j_ss>0 && j_ss<=n
                   Phi(i_ss+1,j_ss+1,1) = ((G(i_ss+1,j_ss+1) + G(i_s+1,j_s+1))/2)*((U(i_s+1,j_s+1) - U(i_ss+1,j_ss+1))/h_1^2);
               end
               i_ss = I_1(4,1); j_ss = I_1(4,2);        % (i+1,j) neighbour.
               if i_ss>0 && i_ss<=m && j_ss>0 && j_ss<=n
                   Phi(i_s+1,j_s+1,2) = ((G(i_s+1,j_s+1) + G(i_ss+1,j_ss+1))/2)*((U(i_ss+1,j_ss+1) - U(i_s+1,j_s+1))/h_2^2);
               end
               i_ss = I_1(5,1); j_ss = I_1(5,2);        % (i-1,j) neighbour.
               if i_ss>0 && i_ss<=m && j_ss>0 && j_ss<=n
                   Phi(i_ss+1,j_ss+1,2) = ((G(i_ss+1,j_ss+1) + G(i_s+1,j_s+1))/2)*((U(i_s+1,j_s+1) - U(i_ss+1,j_ss+1))/h_2^2);
               end
            end
        end

        %% Pixel velocity computation:
        [s1,~] = size(I_1);
        for i = 1:s1
            i_s = I_1(i,1); j_s = I_1(i,2);
            if i_s>0 && i_s<=m && j_s>0 && j_s<=n
                u_dt(i_s,j_s) = ((G(i_s+2,j_s+1) + G(i_s+1,j_s+1))/2)*((U(i_s+2,j_s+1)-U(i_s+1,j_s+1))/(h_1^2)) - ...
               ((G(i_s+1,j_s+1) + G(i_s,j_s+1))/2)*((U(i_s+1,j_s+1)-U(i_s,j_s+1))/(h_1^2)) + ...
               ((G(i_s+1,j_s+2) + G(i_s+1,j_s+1))/2)*((U(i_s+1,j_s+2)-U(i_s+1,j_s+1))/(h_2^2)) - ...
               ((G(i_s+1,j_s+1) + G(i_s+1,j_s))/2)*((U(i_s+1,j_s+1)-U(i_s+1,j_s))/(h_2^2));
                %u_dt(i_s,j_s) = Phi(i_s+1,j_s+1,2) - Phi(i_s,j_s+1,2) + Phi(i_s+1,j_s+1,1) - Phi(i_s+1,j_s,1);
            end
        end

        %% Flow expiry time updates:
        [s2,~] = size(I_2);
        for i = 1:s2
            for j = 1:s2
                i_s = I_2(i,1); j_s = I_2(i,2); i_ss = I_2(j,1); j_ss = I_2(j,2);         
                if i_s>0 && i_s<=m && j_s>0 && j_s<=n && i_ss>0 && i_ss<=m && j_ss>0 && j_ss<=n  
                    if (i_s==i_ss && j_s==j_ss-1) || (i_s==i_ss-1 && j_s==j_ss)
                        tau_s = t_max - time(i_s,j_s);
                        if U(i_s+1,j_s+1)~=U(i_ss+1,j_ss+1) && sign(u_dt(i_s,j_s)-u_dt(i_ss,j_ss))==-sign(U(i_s+1,j_s+1)-U(i_ss+1,j_ss+1))...
                                && abs(U(i_s+1,j_s+1)-U(i_ss+1,j_ss+1))<tau_s*abs(u_dt(i_s,j_s)-u_dt(i_ss,j_ss))
                            tau_s = -(U(i_s+1,j_s+1)-U(i_ss+1,j_ss+1))/(u_dt(i_s,j_s)-u_dt(i_ss,j_ss));
                        end
                        if i_s==i_ss && j_s==j_ss-1
                            data(i_s,j_s,1) = time(i_s,j_s) + tau_s;
                            Q(i_s,j_s,1) = t_max - data(i_s,j_s,1);
                        else
                            data(i_s,j_s,2) = time(i_s,j_s) + tau_s;
                            Q(i_s,j_s,2) = t_max - data(i_s,j_s,2);
                        end
                    end
                end
            end
        end

        %% Pixel expiry time updates:
        [s1,~] = size(I_1);
        for i=1:s1
            i_s = I_1(i,1); j_s = I_1(i,2);
            if i_s>0 && i_s<=m && j_s>0 && j_s<=n
                tau_s = t_max - time(i_s,j_s);
                if u_dt(i_s,j_s)>0
                    [i_tag,j_tag] = FindMaxNeighbour(U(2:m+1,2:n+1),i_s,j_s);
                    if (U(i_s+1,j_s+1) + tau_s*u_dt(i_s,j_s))>U(i_tag+1,j_tag+1)
                        tau_s = (U(i_tag+1,j_tag+1)-U(i_s+1,j_s+1))/u_dt(i_s,j_s);
                    end
                end
                if u_dt(i_s,j_s)<0
                    [i_tag,j_tag] = FindMaxNeighbour(-U(2:m+1,2:n+1),i_s,j_s);
                    if (U(i_s+1,j_s+1) + tau_s*u_dt(i_s,j_s))<U(i_tag+1,j_tag+1)
                        tau_s = (U(i_tag+1,j_tag+1)-U(i_s+1,j_s+1))/u_dt(i_s,j_s);
                    end
                end
                if ( tau_s ~= 0)
                    data(i_s,j_s,3) = time(i_s,j_s) + tau_s;
                end
                Q(i_s,i_s,3) = t_max-data(i_s,j_s,3);
            end
        end

        %fprintf(1, 'cntWhile_nnz = %6d\t', cntWhile_nnz); 
        %fprintf(1, 'cntNumPixUpdate = %3d\t', cntNumPixUpdate); 
        %fprintf(1, 'number of changed values = %3d\t', sum((timePrev(:) - time(:)) ~= 0)); 
        %fprintf(1, 'cntNumOfnoneZeros = %3d\t', nnz(time)); 
        %fprintf(1, 'cntNumNZ = %3d\n', nnz(time-t_max)); 
        %timePrev = time;

    end
end
    
 u = U(2:m+1,2:n+1);
 

                

   
        
        
function u = randFAB_LocallyAdaptiveStepSize(u,t_max,h_1,h_2,g)

[m,n] = size(u);

%% Initialisation:
T = zeros(m,n,2);
T(1:m-1,1:n-1,:) = t_max;   
T(m,1:n-1,1) = t_max;                   % (:,:,1) is the (i,j+1) neighbour.
T(1:m-1,n,2) = t_max;                   % (:,:,2) is the (i+1,j) neighbour.

%% Diffusivity matrix computation:
U = ones(m+2,n+2);
U(2:m+1,2:n+1) = u;
U(1,1) = u(1,1); U(1,2:n+1) = u(1,1:n); U(1,n+2) = u(1,n);
U(2:m+1,1) = u(1:m,1); U(m+2,1) = u(1,m);
U(m+2,2:n+1) = u(m,1:n); U(m+2,n+2) = u(m,n);
U(2:m+1,n+2) = u(1:m,n);

G = ones(m+2,n+2);
for i = 2:m+1
    for j = 2:n+1
        A = max(((U(i+1,j)-U(i,j))/h_1)*((U(i,j)-U(i-1,j))/h_1),0);
        B = max(((U(i,j+1)-U(i,j))/h_2)*((U(i,j)-U(i,j-1))/h_2),0);
        G(i,j) = g(A+B);
    end 
end

%% Asynchronous update:
while nnz(T)~=0
    %% Random selection:
    T_norm = T./max(T(:));
    tau_star = datasample(T(:),1,'Weights',T_norm(:));
    [i,j] = find(T==tau_star);
    [size_1,~] = size(i);
    for i_1=1:size_1
        is = i(i_1,1);
        if j(i_1,1)>n
            js = j(i_1,1) - n;
            is_tag = is + 1;
            js_tag = js; 
            k = 2;
        else
            js = j(i_1,1);
            is_tag = is;
            js_tag = js +1;
            k = 1;
        end

        %% Diffusivity computation:
        g_s = 0.5*(G(is+1,js+1) + G(is_tag+1,js_tag+1));

        %% Flow computation:
       if is == is_tag
           h = h_1;                        % horizontal neighbour.         
       else
           h = h_2;                        % vertical neighbour.   
       end
         u_dt = g_s*(u(is_tag,js_tag)-u(is,js))/h;

       %% Step size determination:
       if g_s > 0 && tau_star > (h^2/(2*g_s))
           tau_star = h^2/(2*g_s);
       end

       [i_star,j_star] = FindMaxNeighbour(u,is,js);
       if g_s < 0 && (is~=i_star || js~=j_star)
           if u(is,js)+tau_star*u_dt > u(i_star,j_star)
               tau_star = (u(i_star,j_star) - u(is,js))/u_dt;
           end
       end

       [i_star,j_star] = FindMaxNeighbour(-u,is,js);
       if g_s < 0 && (is~=i_star || js~=j_star)
           if u(is,js)+tau_star*u_dt < u(i_star,j_star)
               tau_star = (u(i_star,j_star) - u(is,js))/u_dt;
           end
       end

       [i_star,j_star] = FindMaxNeighbour(u,is_tag,js_tag);
       if g_s < 0 && (is_tag~=i_star || js_tag~=j_star)
           if u(is_tag,js_tag)-tau_star*u_dt > u(i_star,j_star)
               tau_star = (u(i_star,j_star) - u(is_tag,js_tag))/(-u_dt);
           end
       end

       [i_star,j_star] = FindMaxNeighbour(-u,is_tag,js_tag);
       if g_s < 0 && (is_tag~=i_star || js_tag~=j_star)
           if u(is_tag,js_tag)-tau_star*u_dt < u(i_star,j_star)
               tau_star = (u(i_star,j_star) - u(is_tag,js_tag))/(-u_dt);
           end
       end

       %% Two pixel flow update:
       u(is,js) = u(is,js) + tau_star*u_dt;
       u(is_tag,js_tag) = u(is_tag,js_tag) - tau_star*u_dt;

       T(is,js,k) = T(is,js,k) - tau_star;      
    end
end
    
    
    
        
       
        
        
    
    
    
    
    
    
    
  
    
    
    
    
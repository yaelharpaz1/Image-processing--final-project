function [u,tau] = FAB_GlobalStepSize(u,t_max,t_min,h_1,h_2,g)

[m,n] = size(u);

%% Diffusivity computation: 
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

%% Flow computation: 
u_dt = ones(m,n);
for i = 2:m+1
    for j = 2:n+1
        u_dt(i-1,j-1) = ((G(i+1,j) + G(i,j))/2)*((U(i+1,j)-U(i,j))/(h_1^2)) - ...
               ((G(i,j) + G(i-1,j))/2)*((U(i,j)-U(i-1,j))/(h_1^2)) + ...
               ((G(i,j+1) + G(i,j))/2)*((U(i,j+1)-U(i,j))/(h_2^2)) - ...
               ((G(i,j) + G(i,j-1))/2)*((U(i,j)-U(i,j-1))/(h_2^2));
    end 
end

%% Step size determination:
tau = t_max;
for i = 1:m
    for j = 1:n  
        [i_tag,j_tag] = FindMaxNeighbour(u,i,j);
        if i~=i_tag || j~=j_tag
           if (u(i,j) + tau*u_dt(i,j)) > (u(i_tag,j_tag) + tau*u_dt(i_tag,j_tag))
               tau_star = -(u(i_tag,j_tag)-u(i,j))/(u_dt(i_tag,j_tag)-u_dt(i,j));
               if tau_star >= t_min
                   tau = tau_star;
               end
           end
        end
        [i_tag,j_tag] = FindMaxNeighbour(-u,i,j);
        if i~=i_tag || j~=j_tag
           if (u(i,j) + tau*u_dt(i,j)) < (u(i_tag,j_tag) + tau*u_dt(i_tag,j_tag))
               tau_star = -(u(i_tag,j_tag)-u(i,j))/(u_dt(i_tag,j_tag)-u_dt(i,j));
               if tau_star >= t_min
                   tau = tau_star;
               end
           end
        end
    end
end

%% Global update:
for i = 1:m
    for j = 1:n  
        u(i,j) = u(i,j) + tau*u_dt(i,j);
    end
end

            
            
        
        
     
        

    


        



        
    

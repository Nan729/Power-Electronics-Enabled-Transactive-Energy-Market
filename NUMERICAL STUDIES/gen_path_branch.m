function [H] = gen_path_branch(branch_info,N)

    M=size(branch_info,1);
    incidence=zeros(N,1);
    branch_end=zeros(N,1);
    for i=1:M
        s=branch_info(i,1);
        e=branch_info(i,2);
        incidence(e)=s;
        branch_end(e)=i;
    end
    search_flag=zeros(N,1);
    search_flag(1,1)=1;
    H=zeros(M,N); 
    H(1,1)=1;
    for k=2:N 
        e=k;
        [H,search_flag] = search_path(H,e,incidence,branch_end,search_flag);
    end
   % H_modified=H(:,2:N);
end


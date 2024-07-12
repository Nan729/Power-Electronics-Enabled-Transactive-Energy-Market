function [H,search_flag] = search_path(H,e,incidence,branch_end,search_flag)

    s=incidence(e);
    if search_flag(s,1)~=1
         [H,search_flag]=search_path(H,s,incidence,branch_end,search_flag);
    end
    H(:,e)=H(:,s);
        % add its own line
    l_relate=branch_end(e);
    H(l_relate,e)=1;  
    search_flag(e,1)=1;
end


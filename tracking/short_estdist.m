function [estdist] = short_estdist (est_dist,nF)

v=0;
for F=1:nF
    if ~isnan(est_dist(F,:))
    v=v+1;
    estdist(v,:)=est_dist(F,:);
    else v=v;
    end
end

end

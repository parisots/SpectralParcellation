function lab = segm(D,faces,Seed,lab,tsh)

lab(Seed) = Seed;
[a,b]=find(faces==Seed);
ConnV=[];
for i=1:numel(a)
ConnV=[ConnV,faces(setdiff(1:3,a(i)),b(i))];
end
ConnV=unique(ConnV(:));
    for j=1:numel(ConnV)
        S=ConnV(j);
        if D(S)<tsh && lab(ConnV(j))==0
            lab(ConnV(j)) = Seed;
            lab = segm(D,faces,ConnV(j),lab,tsh);
        end
    end
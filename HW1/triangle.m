
function n_triangles = triangle(A, isdirected)
    %nodes = size(A,1)
    count_Triangle = 0
    for i = 1:size(A,1)
        for j = 1:size(A,1)
            for k = 1:size(A,1)
                if i~=j & j~=k & k~=i & A(i,j) & A(j,k) & A(k,i)
                    count_Triangle = count_Triangle + 1
                end
            end
        end
    end
    if isdirected
        n_triangles = count_Triangle/3
    else
        n_triangles = count_Triangle/6
    end
end
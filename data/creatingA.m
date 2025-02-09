function A =  creatingA(m)
points = rand(m, 2) * 10;
A = zeros(m);
for i = 1:m
    for j = 1:m
        A(i,j) = norm([points(i,1), points(i,2)]...
                -[points(j,1),points(j,2)]);
    end
end
end
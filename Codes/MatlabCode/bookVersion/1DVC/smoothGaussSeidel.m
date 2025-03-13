function u = smoothGaussSeidel(f,u,m,sweepDirection)

if m == 0
  return
end

nl = length(f);
hl = 1/(nl+1);

for sweep = 1:m
  if strcmp(sweepDirection,'Forward')

    % Update first point.
    u(1) = (hl*f(1)+u(2))/2.0;

    % Update interior points.
    for i = 2:nl-1
      u(i) = (hl*f(i)+u(i-1)+u(i+1))/2.0;
    end

    % Update last point.
    u(nl) = (hl*f(nl)+u(nl-1))/2.0;

  else
    
    % Update last point.
    u(nl) = (hl*f(nl)+u(nl-1))/2.0;

    % Update interior points.
    for i = nl-1:-1:2
      u(i) = (hl*f(i)+u(i-1)+u(i+1))/2.0;
    end

    % Update first point.
    u(1) = (hl*f(1)+u(2))/2.0;

  end
end

end
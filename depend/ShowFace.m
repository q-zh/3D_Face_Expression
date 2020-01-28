function p = ShowFace(input, edge, landmark)
clf;
hold on;

if size(input,1)>100 && size(input,2)>100
    z = input;
    mask = edge;
    [vertex, edge] = find_vertex_edge(z, mask);
    input = vertex;
end


if nargin==1
    if isfield(input, 'normals')
        PlotTarget1 = rmfield(input, 'normals');
    else
        PlotTarget1 = input;
    end
elseif ~isfield(input, 'vertices')
    if size(input,1) == 15285 || size(input,1) ==15282
        input = reshape(input, [3, size(input,1)/3]);
        input = input';
    end
    PlotTarget1.vertices = input;
    PlotTarget1.faces = edge;
    if size(edge,2)~=3
        return;
    end
else
    if isfield(input, 'normals')
        PlotTarget1 = rmfield(input, 'normals');
    else
        PlotTarget1 = input;
    end
end




p = patch(PlotTarget1, 'EdgeColor',  'none');
material([0.5 0.8 0.2])
l = light('Position',[-0.0 0.2 0.9],'Style','infinite');
lighting gouraud
colormap(copper) 
view([0,90]); 
p.FaceColor = [0.5,0.5,0.5];    % set the face colors to be interpolated
axis equal off

hold on
if nargin==3 || (nargin==2 && isfield(input, 'vertices'))
    if nargin == 2
        landmark =edge;
    end
    if size(landmark,2)~=3 &&  size(landmark,1)==1
        plot3(PlotTarget1.vertices(landmark,1), PlotTarget1.vertices(landmark,2), PlotTarget1.vertices(landmark,3), 'r.', 'MarkerSize', 20);
    elseif size(landmark,2)==3
        plot3(landmark(:,1), landmark(:,2), landmark(:,3), 'r.', 'MarkerSize', 20);
    end
        
end
end
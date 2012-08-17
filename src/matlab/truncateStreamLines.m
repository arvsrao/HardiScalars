outdir='/home/arvind/temp/';
tracks = read_mrtrix_tracks([ outdir 'STRtracksFewer.tck' ]);
outdir='/home/arvind/temp/TrackFiles/OK/';
writedir='/home/arvind/MATLAB/StreamFiles/';

for ii=1:100
    if exist( [ outdir 'track_' num2str(ii) '.tck' ] )~=0
       disp( ['processing streamline ' num2str(ii)]);
       line = cell2mat(tracks.data(ii));
       [pp W T] = CutStreamLine( line );

       %qq = ppdiff(pp, 1);
       %dW = ppval(qq, T);
        
       PrintStreamLineFile( W',  [ 'Streamline_bspline_' num2str(ii) '.txt' ] );
       %PrintStreamLineFile( dW', [ 'Streamline_bspline_directions_' num2str(ii) '.txt' ] );
    end
end

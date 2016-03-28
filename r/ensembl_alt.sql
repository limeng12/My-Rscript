select seq_region.name,exon.seq_region_start,exon.seq_region_end
from seq_region inner join exon on seq_region.seq_region_id=exon.seq_region_id
where is_constitutive=0;
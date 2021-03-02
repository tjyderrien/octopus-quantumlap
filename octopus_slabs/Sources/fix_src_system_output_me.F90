diff --git a/src/system/output_me.F90 b/src/system/output_me.F90
index f83dfd6f9..635bba8b2 100644
--- a/src/system/output_me.F90
+++ b/src/system/output_me.F90
@@ -234,7 +234,8 @@ contains
         case(3)
           do ll = 1, this%ks_multipoles
             do mm = -ll, ll
-              write(fname,'(i4)') id
+!              write(fname,'(i4)') id
+               write(fname,'(i7)') id
               write(fname,'(a)') trim(dir)//'/ks_me_multipoles.'//trim(adjustl(fname))
               if (states_are_real(st)) then
                 call doutput_me_ks_multipoles(fname, namespace, st, gr, ll, mm, ik)
@@ -278,7 +279,7 @@ contains
       ASSERT(.not. st%parallel_in_states)
       ! The content of each file should be clear from the header of each file.
       do ik = st%d%kpt%start, st%d%kpt%end
-        write(fname,'(i4)') ik
+        write(fname,'(i7)') ik
         write(fname,'(a)') trim(dir)//'/ks_me_dipole.k'//trim(adjustl(fname))//'_'
           if (states_are_real(st)) then
             call doutput_me_dipole(this, fname, namespace, st, gr, hm, geo, ik)

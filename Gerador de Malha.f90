!  GeradordeMalha.f90 
!
!  FUNCTIONS:
!  GeradordeMalha - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: GeradordeMalha
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program GeradordeMalha
    
    !VARIÁVEIS
    INTEGER NDIME, NNO, NEL, NNOEL, M, Q, C, R, S, W, O !Número de dimensões, Número de nós, Número de elementos, Número de nós por elemento
    REAL*8 L, A, X, Y, PI, da, db, dc, P, B, N, A1, B1, Dab, X1, X2, H, D, E, U, V
    INTEGER, ALLOCATABLE:: INC(:,:), VET(:), IREF(:,:), IREF2(:,:), CONT(:), ELL(:,:), ELT(:,:)
    REAL*8, ALLOCATABLE:: COORD(:,:), REF(:,:), DIST(:)
    
    IL=0
    PI=0.d00
    PI=4.*datan(1.d0)
    P=1
    C=1
    M=0
    S=0 !Parâmetro de contagem do último ponto da fronteira
    
    
    !LEITURA DE DADOS
    OPEN(5,FILE='INPUT.TXT', STATUS='OLD')
    !OPEN(6,FILE='OUTPUTN.TXT', STATUS='UNKNOWN')
    OPEN(6,FILE='OUTPUTN.POS', STATUS='UNKNOWN')
    OPEN(7,FILE='OUTPUTP.TXT', STATUS='UNKNOWN')
    !OPEN(7,FILE='OUTPUTP.POS', STATUS='UNKNOWN')
    READ(5,*) NDIME
    READ(5,*) NNO
    READ(5,*) NEL
    READ(5,*) NNOEL
    ALLOCATE (COORD(NNO,NDIME),INC(NEL,NNOEL), VET(NEL))
    
           
    DO I=1,NNO                                                                                        
       READ(5,*) (COORD(I,J),J=1,NDIME)
    ENDDO       
    
    DO I=1,NEL
        READ(5,*) (INC(I,J),J=1,NNOEL),(VET(I))
    ENDDO
    NR=0
    NR=SUM(VET)
    ALLOCATE (REF(10000,2))
    ALLOCATE (IREF2(10000,3))
    ALLOCATE (IREF(10000,3))
    ALLOCATE (CONT(5))
    ALLOCATE (ELL(10000,3))
    ALLOCATE (ELT(10000,3))
    
    ! CORPO DO PROGRAMA
    
    !***********************************************************************************************************************************************
    !***********************************************************************************************************************************************
                                                                    !OUTPUT NAESY
    !***********************************************************************************************************************************************
    !***********************************************************************************************************************************************                                                                
                                                              !Matriz de Contorno Refinada
    DO I=1,NEL
        A = (2.d+00)/(VET(I))
        n=-1d+00
        DO B=1,((VET(I)))        
            da=0.5d+00*n*(n-1.d+00) !Funções de Forma
            db=(1.d+00-n)*(1.d+00+n) !Funções de Forma
            dc=0.5d+00*n*(n+1.d+00) !Funções de Forma
            IL=IL+1
            x=da*(COORD(INC(I,1),1))+db*(COORD(INC(I,2),1))+dc*(COORD(INC(I,3),1))               ! Verificar
            IF (B<=NR) THEN
                REF(IL,1) = da*(COORD(INC(I,1),1))+db*(COORD(INC(I,2),1))+dc*(COORD(INC(I,3),1))    !COORDENADAS X
                REF(IL,2) = da*(COORD(INC(I,1),2))+db*(COORD(INC(I,2),2))+dc*(COORD(INC(I,3),2))     !COORDENADAS Y
                n=-1+(A*B)
            ELSE 
            ENDIF
            
        ENDDO
    ENDDO  
    
    
    !Repetição do ponto 1 para o ponto NR+1
    S3=1
    S1=NR
    S2=S1+1
    S4=S2+1
    !REF(S2,1)=COORD(INC(1,1),1)
    !REF(S2,2)=COORD(INC(1,1),2)
    
    !***********************************************************************************************************************************************
                                                                    !Matriz de Incidência de Contorno Refinada
    B=2
    DO I=1, NR
        B=B-1
        DO J=1,NDIME
            IF (I==NR) THEN
                IF (J==2) THEN
                    B=1
                ENDIF
            ENDIF
            IREF(I,J)=B
            B=B+1   
        ENDDO
    ENDDO       
    !***********************************************************************************************************************************************
                                                                    !ARQUIVO .TXT NAESY
    WRITE (6,'(A)') '%HEADER'
    WRITE (6,'(A)') 'NEUTRAL FILE CREATED BY FEMOOP PROGRAM'
    WRITE (6,*)
    
    WRITE (6,'(A)') '%NODE'
    WRITE (6,*) NR
    WRITE (6,*)
    
    WRITE (6,'(A)') '%NODE.COORD'
    WRITE (6,*) NR
        DO I=1,NR
            WRITE(6,'(I6,2X,f18.7,f25.7,f25.7)') I, (REF(I,J), J=1,NDIME), 0
        ENDDO  
    WRITE(6,*)
    WRITE(6,*)' %INTEGRATION.QUADRATURE '
    WRITE(6,*) ' 2 '
    WRITE(6,*)
    WRITE(6,*)' %INTEGRATION.QUADRATURE.QUADRILATERAL'
    WRITE(6,*) '1 '
    WRITE(6,*) '1    1   1   1  1 ' 
    WRITE(6,*)
    WRITE(6,*)'%ELEMENT'
    WRITE(6,*) NR
    
    WRITE(6,*)
    WRITE(6,*) '%ELEMENT.T3'
    WRITE(6,*) NR
        DO  I=1,NR                        
            write(6,'(i6,2x,i3,2x,i3,2x,i3,2x,i4,2x,i4,2x,i4)') I, 1, 1, 1, (IREF(I,J),J=1,NDIME), (IREF(I,2))
        ENDDO 
    WRITE (6,'(A)')'%END'
    
    !Criação da matriz de elementos lineares ELL
    DO I=1,NR
        ELL(I,1)=IREF(I,1)
        ELL(I,2)=IREF(I,2)
        ELL(I,3)=0
    ENDDO
    
    
    !***********************************************************************************************************************************************
    !***********************************************************************************************************************************************
                                                                    !OUTPUT POS3D
    !***********************************************************************************************************************************************
    !***********************************************************************************************************************************************                                                               
                                                            !Malha Interna De Fronteira 
    P=0
    I=1  !ELEMENTO lINEAR
    J=0  !ELEMENTO TRIANGULAR
    W=NR !NÚMERO DE ELEMENTOS 
    R=NR
    Q=0  !NÚMERO DE PONTOS
    O=0
    S1=S1+1
    C=1
    F=0
    B=0
    
    DO WHILE (P==0)
        
            DO
                C=C+1
                G=ELL(C,3)
                IF (C==R+1) EXIT
                IF (G==0) EXIT
            ENDDO
        
            G=ELL(C,1)
            ELL(C,3)=1
            
            DO
                B=B+1
                N=ELL(B,2)
                T=ELL(B,3)
                IF (B==R+1) EXIT
                IF (T==1) CYCLE
                IF (N==G) EXIT
                IF (N/=G) CYCLE
                
            ENDDO
        Z=0
        IF (B==R+1) THEN
            Z=1
            B=0
            ELL(C,3)=0
        ENDIF
        IF (Z==0) THEN 
         A1=SQRT(((REF(ELL(C,2),1)-REF(ELL(C,1),1))**2)+((REF(ELL(C,2),2)-REF(ELL(C,1),2))**2))   !Norma do elemento I+1
         A2=SQRT(((REF(ELL(B,2),1)-REF(ELL(B,1),1))**2)+((REF(ELL(B,2),2)-REF(ELL(B,1),2))**2))   !Norma do elemento B
         
         V1=REF(ELL(C,2),1)-REF(ELL(C,1),1)   !componetne vetorial x do elemento I+1 (X)
         V2=REF(ELL(C,2),2)-REF(ELL(C,1),2)   !componente vetorial y do elemento I+1 (Y)
         
         V3=REF(ELL(B,1),1)-REF(ELL(B,2),1)   !componetne vetorial x do elemento B (X)
         V4=REF(ELL(B,1),2)-REF(ELL(B,2),2)   !componente vetorial y do elemento B (Y) 
         
         A=ACOS((V1*V3+V2*V4)/(A1*A2)) ! Ângulo entre dois vetores
         
         !Verificação de vértice
         A3=0.5*(REF(ELL(B,1),1)*REF(ELL(B,2),2)-REF(ELL(B,1),2)*REF(ELL(B,2),1)+REF(ELL(C,1),1)*REF(ELL(C,2),2)-REF(ELL(C,1),2)*REF(ELL(C,2),1)) 
         
         IF (A3<=0) THEN
             Y=Y+1
         ENDIF
         
         !IF (A3<0) THEN
            ! A=PI-2*A       !Correção do ângulo interno
         !ENDIF
         
         IF (A<=PI/2) THEN
             
            J=J+1
            IREF2(J,1)=ELL(B,1)
            IREF2(J,2)=ELL(B,2)
            IREF2(J,3)=ELL(C,2)
            
            R=R+1
            ELL(R,1)=ELL(B,1)
            ELL(R,2)=ELL(C,2)
            ELL(R,3)=0
            
            
            
            ELL(B,3)=1
            ELL(C,3)=1
            
            O=O+1

         ELSEIF (PI/2<A<=3*PI/2) THEN
             H1=SQRT(3.0)*A2/2   !ALTURA DO TRIÂNGULO
             H2=SQRT(3.0)*A1/2
             H3=0.86
             
             
             A=A/2 !ÂNGULO INTERNO
             B1=ACOS(V1/A1)  !ÂNGULO EM RELAÇÃO AO VETOR UNITÁRIO X
             
             C1=REF(ELL(C,2),2)-REF(ELL(C,1),2)
             IF (C1<0) THEN
                 B1=B1*(-1)
             ENDIF
             
             A=B1+A !SOMA DOS ÂNGULOS
             REF(S1,1)=REF(ELL(C,1),1)+H3*COS(A)  !SOMA DAS COORDENADAS
             REF(S1,2)=REF(ELL(C,1),2)+H3*SIN(A)  !SOMA DAS COORDENADAS
             
             Q=Q+1
             J=J+1  
             IREF2(J,1)=ELL(B,1)
             IREF2(J,2)=ELL(B,2)
             IREF2(J,3)=S1
             
             J=J+1
             IREF2(J,1)=ELL(C,1)
             IREF2(J,2)=ELL(C,2)
             IREF2(J,3)=S1
             
             R=R+1
             ELL(R,1)=ELL(B,1)
             ELL(R,2)=S1
             ELL(R,3)=0
             
             R=R+1
             ELL(R,1)=S1
             ELL(R,2)=ELL(C,2)
             ELL(R,3)=0
             
             ELL(B,3)=1
             ELL(C,3)=1
             
             S1=S1+1
             O=O+2
         ELSEIF (3*PI/2<A<=2*PI) THEN
            H2=SQRT(3.0)*A2/2
            
            B1=ACOS(V3/A2) !Ângulo em relação ao vetor unitário X
            
            C1=REF(ELL(B,2),2)-REF(ELL(B,2),2)
            IF (C1<0) THEN
                 B1=B1*(-1)
            ENDIF
            
            A=B1+PI/2 !Ângulo total
            REF(S1,1)=(0.5*REF(ELL(B,1),1)+0.5*REF(ELL(B,2),1))+H2*COS(A)
            REF(S1,2)=(0.5*REF(ELL(B,1),2)+0.5*REF(ELL(B,2),2))+H2*SIN(A)
            
            Q=Q+1
            J=J+1
            IREF2(J,1)=ELL(B,1)
            IREF2(J,2)=ELL(B,2)
            IREF2(J,3)=S1
            
            R=R+1
            ELL(R,1)=ELL(B,1)
            ELL(R,2)=S1
            ELL(R,3)=0
            
            R=R+1
            ELL(R,1)=S1
            ELL(R,2)=ELL(B,2)
            ELL(R,3)=0
            
            ELL(B,3)=0
            ELL(C,3)=0
            
            S1=S1+1
            O=O+1
            
         ENDIF
         
         
         C=1
         B=0
         I=I+1
         IF (I>=50) THEN   
         P=P+1
         ENDIF
        ENDIF
    ENDDO
    

    
    !***********************************************************************************************************************************************
                                                                    !ARQUIVO .TXT POS3D
                                                                    
    WRITE (7,'(A)') '%HEADER'
    WRITE (7,'(A)') 'NEUTRAL FILE CREATED BY FEMOOP PROGRAM'
    WRITE (7,*)
    
    WRITE (7,'(A)') '%NODE'
    WRITE (7,*) NR+Q
    WRITE (7,*)
    
    WRITE (7,'(A)') '%NODE.COORD'
    WRITE (7,*) NR+Q
        DO I=1,NR+Q
            WRITE(7,'(I6,2X,f18.7,f25.7,f25.7)') I, (REF(I,J), J=1,NDIME), 0
        ENDDO  
    WRITE(7,*)
    WRITE(7,*)' %INTEGRATION.QUADRATURE '
    WRITE(7,*) ' 2 '
    WRITE(7,*)
    WRITE(7,*)' %INTEGRATION.QUADRATURE.QUADRILATERAL'
    WRITE(7,*) '1 '
    WRITE(7,*) '1    1   1   1  1 ' 
    WRITE(7,*)
    WRITE(7,*)'%ELEMENT'
    WRITE(7,*) O
    
    WRITE(7,*)
    WRITE(7,*) '%ELEMENT.T3'
    WRITE(7,*) O
        DO  I=1,O                  
            write(7,'(i6,2x,i3,2x,i3,2x,i3,2x,i4,2x,i4,2x,i4)') I, 1, 1, 1, (IREF2(I,J),J=1,3)
        ENDDO 
    WRITE (7,'(A)')'%END'                                                              
                                                                    
                                                                    
                                                                    
    !***********************************************************************************************************************************************
    !***********************************************************************************************************************************************    
    CONTAINS 
   
    SUBROUTINE AJUSTE
         N=ACOS(((X1*1)+(Y1*0))/(A*1))                   !Ângulo entre o vetor e o vetor unitário EMBUTIDO
         IF (REF(IREF(I,2),2)-REF(IREF(I,1),2) < 0) THEN    !
            N=N*(-1)                                        !
            ELSE                                            ! Ajusta de direção
            N=N                                             !
        ENDIF                                               !
    END SUBROUTINE AJUSTE
    
    SUBROUTINE AJUSTE1
         N=ACOS(((X1*1)+(Y1*0))/(A*1))                   !Ângulo entre o vetor e o vetor unitário EMBUTIDO
         IF (REF(IREF2(W+I,2),2)-REF(IREF2(W+I,1),2) < 0) THEN    !
            N=N*(-1)                                        !
            ELSE                                            ! Ajusta de direção
            N=N                                             !
        ENDIF                                               !
    END SUBROUTINE AJUSTE1
    
    SUBROUTINE AREA
        P=0
        B=0
        DO I=1,NR

           M=I+1
           IF (M .GE. NR+1) THEN
               M=1
           ENDIF 
           P=0
           P=((REF(NR+I,1)*REF(M,2))-(REF(NR+I,2)*REF(M,1)))
           B=B+P
        ENDDO
        P=0
        P=B/2
    END SUBROUTINE
       
    SUBROUTINE ALTURA
        IF ((D/H)>=0.8 .AND. (D/H)<=1.2) THEN
            H=H
        ELSE
            H=D
        ENDIF
    END SUBROUTINE
   
    SUBROUTINE DISTA
        DO I=47,(NR)
         DO J=I+1,NR
          A=SQRT((((REF(IREF2(J,3),1))-(REF(IREF2(I,3),1)))**2)+(((REF(IREF2(J,3),2))-(REF(IREF2(I,3),2)))**2))
         ENDDO
       ENDDO   
    END SUBROUTINE
    
    SUBROUTINE PIN
        !Ponto no interior da fronteira
        DO O=1,NR
            X3=D-REF(IREF(O,1),1)    ! Vetor u
            Y3=E-REF(IREF(O,2),2)    !   
            
            X4=REF(IREF(O,2),1)-REF(IREF(O,1),1)    ! Vetor v
            Y4=REF(IREF(O,2),2)-REF(IREF(O,1),2)    !
        
            Z= ((-X3*Y4)+(X4*Y3)) !Coordenada Z do vetor ortogonal
            
            IF (Z<0) THEN
                F=1
            ENDIF  
        ENDDO
    END SUBROUTINE
    
   
                                                                      
   end program GeradordeMalha


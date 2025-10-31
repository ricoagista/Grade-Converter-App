import java.util.Scanner;

public class tugas1 {
    
    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        System.out.print("Masukkan jumlah angka: ");
        int n = scanner.nextInt();
        
        int[] angka = new int[n];
        System.out.println("Masukkan angka-angka:");
        for (int i = 0; i < n; i++) {
            angka[i] = scanner.nextInt();
        }
        
        int total = 0;
        for (int i = 1; i < n; i++) { // Perhatikan kesalahan pada index perulangan
            total += angka[i];
        }
        
        double rataRata = (double) total / n;
        System.out.println("Rata-rata adalah: " + rataRata);
        
        scanner.close();
    }
}

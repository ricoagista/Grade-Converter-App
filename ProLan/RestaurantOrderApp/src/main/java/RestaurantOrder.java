import java.util.Scanner;

public class RestaurantOrder {
    // Atribut
    private String customerName;
    private String menuName;
    private int quantity;

    private double pricePerItem;

    // Konstruktor
    public RestaurantOrder(String customerName, String menuName, int quantity, double pricePerItem) {
        this.customerName = customerName;
        this.menuName = menuName;
        this.quantity = quantity;
        this.pricePerItem = pricePerItem;
    }

    // Method menghitung total harga
    public double calculateTotal() {
        // Custom Live Template bisa dibuat untuk rumus ini, misalnya dengan abbreviation: calcTotal
        return quantity * pricePerItem;
    
    }

    // Method menampilkan nota
    public void printReceipt() {
        System.out.println("\n====== NOTA PEMESANAN ======");
        System.out.println("Nama Pelanggan : " + customerName);
        System.out.println("Menu           : " + menuName);
        System.out.println("Jumlah Pesanan : " + quantity);
        System.out.println("Harga Satuan   : Rp" + pricePerItem);
        System.out.println("----------------------------");
        System.out.println("Total Harga    : Rp" + calculateTotal());
        System.out.println("============================");
    }

    // Main Program
    public static void main(String[] args) {
        Scanner input = new Scanner(System.in);

        System.out.println("=== Aplikasi Nota Pemesanan Restoran ===");
        System.out.print("Masukkan nama pelanggan: ");
        String name = input.nextLine();

        System.out.print("Masukkan nama menu: ");
        String menu = input.nextLine();

        System.out.print("Masukkan jumlah pesanan: ");
        int qty = input.nextInt();

        System.out.print("Masukkan harga per item: ");
        double price = input.nextDouble();

        RestaurantOrder order = new RestaurantOrder(name, menu, qty, price);
        order.printReceipt();
    }
}
